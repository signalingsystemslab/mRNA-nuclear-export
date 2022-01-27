#!/bin/python3
import argparse, re
          
def main():
    parser = argparse.ArgumentParser(description='Separate gtf into two base on gene_ids.')
    parser.add_argument('--gtf',  help='Path to original gtf file.', type=str, nargs=1, required=True)
    parser.add_argument('--out_gtf', help='Path of the ouput gtf file for genes matching the gene list (default: GTF_in_GENE_LIST.gtf).', type=str, nargs=1)
    parser.add_argument('--gene_list', help='Path of file containing gene ids, gene names, transcripts to keep, start and end of region', type=str, nargs=1, required=True)

    args = parser.parse_args()

    gtf=args.gtf[0]
    gene_list_path = args.gene_list[0]

    if not args.out_gtf:
        in_gtf=gtf.split('.gtf')[0]+'_in_' + gene_list_path.split('.txt')[0] + '.gtf'    
    else:
        in_gtf=args.out_gtf[0]
    
    g = open(args.gene_list[0], 'r')
    
    gene_dict = {'gene_id':[], 'gene_name':[], 'transcripts':[],'transcripts_TSS':[], 'transcripts_TES':[], 'region':[]}
    
    gene_dict_idx=-1
    for gl in g:
        if (gl.startswith('#')):
            continue
        gl_line = gl.strip('\n').split('\t')
        
        tx = []
        if gl_line[0] == '': 
            gene_dict['transcripts'][gene_dict_idx].append(gl_line[8])
            gene_dict['transcripts_TSS'][gene_dict_idx].append(gl_line[6])
            gene_dict['transcripts_TES'][gene_dict_idx].append(gl_line[7])
        else:
            gene_dict_idx += 1
            gene_dict['gene_id'].append(gl_line[0])
            gene_dict['gene_name'].append(gl_line[1])
            gene_dict['transcripts'].append([gl_line[8]])
            gene_dict['transcripts_TSS'].append([gl_line[6]])
            gene_dict['transcripts_TES'].append([gl_line[7]])
            gene_dict['region'].append([gl_line[4], gl_line[5]])
    
    f = open(gtf, 'r')
    f_in = open(in_gtf, 'a')
    
    
    res_gene_dict = {'chr':[], 'source':[], 'feature':[], 'start':[], 'end':[], 'score':[], 'strand':[], 'frame':[], 'attributes':[], 'gene_id':[], 'gene_symbol':[], 'transcript_dict':[]}
    for line in f:
        if (line.startswith('#')):
            continue
        
        line_elem=line.strip('\n').split('\t')    
        
        if line_elem[2] == 'gene':
            gene_id = re.sub('"', '', [s for s in line_elem[8].split('; ') if 'gene_id'in s][0].split(' ')[1])
            # print(gene_id)
            
            if gene_id in gene_dict['gene_id']:
                gene_idx = [i for i in range(len(gene_dict['gene_id'])) if gene_dict['gene_id'][i] == gene_id][0]
            else:
                gene_idx = []

        if line_elem[2] == 'transcript':
            tx_id = re.sub('"', '', [s for s in line_elem[8].split('; ') if 'transcript_id'in s][0].split(' ')[1])
            strand = line_elem[6]
            if strand == "+":
                tss = min(line_elem[3],line_elem[4])
                tes = max(line_elem[3],line_elem[4])
            elif strand == "-":
                tss = max(line_elem[3],line_elem[4])
                tes = min(line_elem[3],line_elem[4])
                
        if gene_id in gene_dict['gene_id'] and line_elem[2] == 'gene':
            if res_gene_dict['gene_id']:
                write_previous_gene(res_gene_dict, f_in)
                
            res_gene_dict = {'chr':[], 'source':[], 'feature':[], 'start':[], 'end':[], 'score':[], 'strand':[], 'frame':[], 'attributes':[], 'gene_id':[], 'gene_symbol':[], 'transcript_dict':[]}

            res_gene_dict['chr'] = line_elem[0]
            res_gene_dict['source'] = line_elem[1]
            res_gene_dict['feature'] = line_elem[2]
            res_gene_dict['start'] =line_elem[3]
            res_gene_dict['end'] = line_elem[4]
            res_gene_dict['score'] = line_elem[5]
            res_gene_dict['strand'] = line_elem[6]
            res_gene_dict['frame'] = line_elem[7]
            res_gene_dict['attributes'] = line_elem[8]
            res_gene_dict['gene_id'] =  gene_id
            res_gene_dict['gene_symbol'] =  gene_dict['gene_name'][gene_idx]
            res_gene_dict['transcript_dict'] = {'transcript_id': [], 'transcript_info': []}

        
        if gene_id in gene_dict['gene_id'] and tx_id in gene_dict['transcripts'][gene_idx]:
            tx_idx = [i for i in range(len(gene_dict['transcripts'][gene_idx])) if gene_dict['transcripts'][gene_idx][i] == tx_id][0]
            if line_elem[2] == 'transcript':
                res_gene_dict['transcript_dict']['transcript_id'].append(tx_id)
                res_gene_dict['transcript_dict']['transcript_info'].append({'chr':[], 'source':[], 'feature':[], 'start':[], 'end':[], 'score':[], 'strand':[], 'frame':[], 'attributes':[]})

            tx_res_idx = [i for i in range(len(res_gene_dict['transcript_dict']['transcript_id'])) if res_gene_dict['transcript_dict']['transcript_id'][i] == tx_id][0]

            if line_elem[2] == 'exon':
                if line_elem[3] == tss:
                    if not gene_dict['transcripts_TSS'][gene_idx][tx_idx] == '?':
                        line_elem[3] = gene_dict['transcripts_TSS'][gene_idx][tx_idx]
                elif line_elem[3] == tes:
                    if not gene_dict['transcripts_TES'][gene_idx][tx_idx] == '?':
                        line_elem[3] = gene_dict['transcripts_TES'][gene_idx][tx_idx]
                elif line_elem[4] == tss:
                    if not gene_dict['transcripts_TSS'][gene_idx][tx_idx] == '?':
                        line_elem[4] = gene_dict['transcripts_TSS'][gene_idx][tx_idx]
                elif line_elem[4] == tes:
                    if not gene_dict['transcripts_TES'][gene_idx][tx_idx] == '?':
                        line_elem[4] = gene_dict['transcripts_TES'][gene_idx][tx_idx]
                
                to_write = False
                if strand == '+':
                    if not (int(line_elem[4]) <= int(gene_dict['region'][gene_idx][0]) or int(line_elem[3]) >= int(gene_dict['region'][gene_idx][1])):
                        if int(line_elem[3]) < int(gene_dict['region'][gene_idx][0]):
                            line_elem[3] = gene_dict['region'][gene_idx][0]
                        if int(line_elem[4]) > int(gene_dict['region'][gene_idx][1]):
                            line_elem[4] = gene_dict['region'][gene_idx][1]
                        line_elem[8] = 'gene_id '+gene_id+'; transcript_id '+tx_id+';'
                        to_write = True
                        # f_in.write('\t'.join(line_elem)+'\n')
                        
                elif strand == '-': 
                    if not (int(line_elem[3]) >= int(gene_dict['region'][gene_idx][0]) or int(line_elem[4]) <= int(gene_dict['region'][gene_idx][1])):
                        if int(line_elem[4]) > int(gene_dict['region'][gene_idx][0]):
                            line_elem[4] = gene_dict['region'][gene_idx][0]
                        if int(line_elem[3]) < int(gene_dict['region'][gene_idx][1]):
                            line_elem[3] = gene_dict['region'][gene_idx][1]
                        line_elem[8] = 'gene_id '+gene_id+'; transcript_id '+tx_id+';'
                        to_write = True
                        # f_in.write('\t'.join(line_elem)+'\n')
                
                if to_write:
                    res_gene_dict['transcript_dict']['transcript_info'][tx_res_idx]['chr'].append(line_elem[0])
                    res_gene_dict['transcript_dict']['transcript_info'][tx_res_idx]['source'].append(line_elem[1])
                    res_gene_dict['transcript_dict']['transcript_info'][tx_res_idx]['feature'].append(line_elem[2])
                    res_gene_dict['transcript_dict']['transcript_info'][tx_res_idx]['start'].append(line_elem[3])
                    res_gene_dict['transcript_dict']['transcript_info'][tx_res_idx]['end'].append(line_elem[4])
                    res_gene_dict['transcript_dict']['transcript_info'][tx_res_idx]['score'].append(line_elem[5])
                    res_gene_dict['transcript_dict']['transcript_info'][tx_res_idx]['strand'].append(line_elem[6])
                    res_gene_dict['transcript_dict']['transcript_info'][tx_res_idx]['frame'].append(line_elem[7])
                    res_gene_dict['transcript_dict']['transcript_info'][tx_res_idx]['attributes'].append(line_elem[8])
    
    # write last gene
    write_previous_gene(res_gene_dict, f_in)
    
    f.close()
    f_in.close()

def write_previous_gene(gene_dict, f):
    # if multiple transcripts -> compare and merge
    if len(gene_dict['transcript_dict']['transcript_id']) > 1 :
        # print(gene_dict)

        exons_start_transcript = []
        exons_end_transcript = []
        
        exons_start_transcript1 = gene_dict['transcript_dict']['transcript_info'][0]['start']
        exons_end_transcript1 = gene_dict['transcript_dict']['transcript_info'][0]['end']
        exons_start_transcript = []
        exons_end_transcript = []

        for i in range(1,len(gene_dict['transcript_dict']['transcript_id'])):
            # print(exons_start_transcript)
            # print(exons_end_transcript)
            
            exons_start_transcript2 = gene_dict['transcript_dict']['transcript_info'][i]['start']
            exons_end_transcript2 = gene_dict['transcript_dict']['transcript_info'][i]['end']
            
            if i > 1: 
                exons_start_transcript1 = exons_start_transcript
                exons_end_transcript1 = exons_end_transcript
                exons_start_transcript = []
                exons_end_transcript = []

            # check exons t2 vs exons t1 put it in t
            for j in range(len(exons_start_transcript2)):
                for l in range(len(exons_start_transcript1)):
                    # overlapping
                    if int(exons_start_transcript2[j]) <= int(exons_end_transcript1[l]) and int(exons_end_transcript2[j]) >= int(exons_start_transcript1[l]):
                        exons_start_transcript.append(max(int(exons_start_transcript2[j]), int(exons_start_transcript1[l])))
                        exons_end_transcript.append(min(int(exons_end_transcript2[j]), int(exons_end_transcript1[l])))
      
        chrom = gene_dict['chr']
        source = gene_dict['source']
        score = gene_dict['score']
        strand = gene_dict['strand']
        frame = gene_dict['frame']
        gene_id = gene_dict['gene_id']
        
        line_to_write = chrom + '\t' + source +'\t' + gene_dict['feature'] + '\t' + gene_dict['start'] + '\t'+gene_dict['end'] + '\t' + score +'\t' + strand + '\t' + frame + '\t' + gene_dict['attributes']+'\n'
        
        # print(exons_start_transcript)
        # print(exons_end_transcript)
        
        for i in range(len(exons_end_transcript)):
            line_to_write += chrom + '\t' + source +'\t' + 'exon' + '\t' + str(exons_start_transcript[i]) + '\t' + str(exons_end_transcript[i]) + '\t' + score +'\t' + strand + '\t' + frame + '\tgene_id ' + gene_id + '; transcript_id ' + ', '.join(gene_dict['transcript_dict']['transcript_id']) +';\n'
        # print(line_to_write)
        
    elif len(gene_dict['transcript_dict']['transcript_id']) == 1:
        exons_start_transcript = gene_dict['transcript_dict']['transcript_info'][0]['start']
        exons_end_transcript = gene_dict['transcript_dict']['transcript_info'][0]['end']
        
        line_to_write = gene_dict['chr'] + '\t' + gene_dict['source'] +'\t' + gene_dict['feature'] + '\t' + gene_dict['start'] + '\t' + gene_dict['end'] + '\t' + gene_dict['score'] +'\t' + gene_dict['strand'] + '\t' + gene_dict['frame'] + '\t' + gene_dict['attributes']+'\n'
        
        tx_dict = gene_dict['transcript_dict']['transcript_info'][0]
        for i in range(len(tx_dict['start'])):
            line_to_write += tx_dict['chr'][i] + '\t' + tx_dict['source'][i] +'\t' + 'exon' + '\t' + tx_dict['start'][i] + '\t' + tx_dict['end'][i] + '\t' + tx_dict['score'][i] +'\t' + tx_dict['strand'][i] + '\t' + tx_dict['frame'][i] + '\t' +  tx_dict['attributes'][i] + '\n'        
    
    # no exons (i.e. too different ends empty 5kb regions)
    else:
        exons_end_transcript == []
    
    # calculate overlapping regions
    overlap = 0
    for i in range(len(exons_end_transcript)):
        overlap += int(exons_end_transcript[i]) - int(exons_start_transcript[i]) + 1
    
    # print gene_symbol and overlap
    print(gene_dict['gene_symbol']+': '+str(overlap)+'bp')
    # write gene and exons if overlapping region > 500
    if overlap >= 500:
        f.write(line_to_write)
        

main()

