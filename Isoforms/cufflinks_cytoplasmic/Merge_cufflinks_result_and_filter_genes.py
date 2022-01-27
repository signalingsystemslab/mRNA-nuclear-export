#!/bin/pyhton3 
import os
import glob
import numpy as np

def main():
    # parser = argparse.ArgumentParser(description='Merge result of cufflinks for table')
    # parser.add_argument('--cufflinks_results',  help='Path to cufflinks results.', type=str,required=True)
    # parser.add_argument('-o', '--out', help='Path of the ouput txt file (default: merged_cufflinks.txt).', type=str, nargs=1)
    # parser.add_argument('--gene_list', help='Path of file containing gene ids for which to get cufflinks results', type=str, nargs=1)

    # args = parser.parse_args()

    # files=args.gtf[0]
    # files = glob.glob('cytoplasmic.LPA.*/*/isoforms.fpkm_tracking', recursive=True)
    files = glob.glob('cytoplasmic.Naive.*/*/isoforms.fpkm_tracking', recursive=True)

    # if not args.out:
    #     out='merged_cufflinks.txt' 
    # else:
    #     out=args.out[0]

    # out = 'cufflink_LPA_merged.txt'
    out = 'cufflink_Naive_merged.txt'

    # if not args.gene_list:
    #     gene_list=[]
    # else:
    
    gene_list = 'gene_list_with_strand.txt'
    g=open(gene_list, 'r')
    gene_list=[gene.rstrip('\n').split('\t') for gene in g.readlines()]
    gene_ids = [g[0] for g in gene_list] 
    gene_strand = [g[1] for g in gene_list] 

    gene_dict = {'transcript_id':[], 'gene_id':[], 'gene_symbol':[], 'locus':[], 'TSS':[],'TES':[], 'length':[], 'Samples':[], 'FPKM':[]}

    for fpath in sorted(files): 
        print(fpath)
        f=open(fpath,'r')
        for line in f:
            line_elem=line.split('\t')
            if line_elem[3] in gene_ids:
                if line_elem[0] in gene_dict['transcript_id']:
                    idx = [i for i in range(len(gene_dict['transcript_id'])) if gene_dict['transcript_id'][i] == line_elem[0]][0]
                    gene_dict['Samples'][idx] += [fpath.split('/')[len(fpath.split('/'))-2]]
                    gene_dict['FPKM'][idx] += [line_elem[9]]
                else:
                    gene_dict['transcript_id'] += [line_elem[0]]
                    gene_dict['gene_id'] += [line_elem[3]]
                    gene_dict['gene_symbol'] += [line_elem[4]]
                    gene_dict['locus'] += [line_elem[6]]
                    strand = [gene_strand[i] for i in range(len(gene_ids)) if gene_ids[i] == line_elem[3]]
                    if strand[0] == '+':
                        gene_dict['TSS'] += [line_elem[6].split(':')[1].split('-')[0]]
                        gene_dict['TES'] += [line_elem[6].split(':')[1].split('-')[1]]
                    else:
                        gene_dict['TSS'] += [line_elem[6].split(':')[1].split('-')[1]]
                        gene_dict['TES'] += [line_elem[6].split(':')[1].split('-')[0]]
                    gene_dict['length'] += [line_elem[7]]
                    gene_dict['Samples'] += [[fpath.split('/')[len(fpath.split('/'))-2]]]
                    gene_dict['FPKM'] += [[line_elem[9]]]
            else:
                continue

    # write and calculate % to export
    fout=open(out, 'w')
    fout.write('transcript_id\tgene_id\tgene_symbol\tlocaus\tTSS\tTES\tlength\tNb_samples_over_10%\tNb_samples_over_25%\tNb_samples_over_50%\tNb_samples_over_75%\t'+'\t'.join(gene_dict['Samples'][0])+'\n')

    for gene in gene_ids:
        ids = [i for i in range(len(gene_dict['gene_id'])) if gene_dict['gene_id'][i] == gene]
        sample_total = [sum([float(gene_dict['FPKM'][i][si]) for i in ids]) for si in range(len(gene_dict['Samples'][ids[0]]))]
        main1_fpkm = [sorted([float(gene_dict['FPKM'][i][si]) for i in ids], reverse=True)[0] for si in range(len(gene_dict['Samples'][ids[0]]))]
        main1_pct = [ main1_fpkm[i] / sample_total[i] if sample_total[i] > 0 else float('nan') for i in range(len(gene_dict['Samples'][ids[0]])) ]
        
        if len(ids) > 1:
            main2_fpkm =  [sorted([float(gene_dict['FPKM'][i][si]) for i in ids], reverse=True)[1] for si in range(len(gene_dict['Samples'][ids[0]]))]
            main2_over_main1 = [ main2_fpkm[i] / main1_fpkm[i] if main1_fpkm[i] > 0 else float('nan') for i in range(len(gene_dict['Samples'][ids[0]]))]
        
        nb_sample_above_10 = [sum([1 for si in range(len(gene_dict['Samples'][ids[0]])) if float(gene_dict['FPKM'][i][si]) > 0.1 * sample_total[si] and sample_total[si] >= 1]) for i in ids]
        nb_sample_above_25 = [sum([1 for si in range(len(gene_dict['Samples'][ids[0]])) if float(gene_dict['FPKM'][i][si]) > 0.25 * sample_total[si] and sample_total[si] >= 1]) for i in ids]
        nb_sample_above_50 = [sum([1 for si in range(len(gene_dict['Samples'][ids[0]])) if float(gene_dict['FPKM'][i][si]) > 0.5 * sample_total[si] and sample_total[si] >= 1]) for i in ids]
        nb_sample_above_75 = [sum([1 for si in range(len(gene_dict['Samples'][ids[0]])) if float(gene_dict['FPKM'][i][si]) > 0.75 * sample_total[si] and sample_total[si] >= 1]) for i in ids]
        total_above_1_fpkm = sum([1 for si in range(len(gene_dict['Samples'][ids[0]])) if sample_total[si] >= 1])
        
        for i in ids:
            i2 = [j for j in range(len(ids)) if ids[j] == i][0]
            line_to_write = [gene_dict['gene_id'][i]]
            line_to_write += [gene_dict['transcript_id'][i]]
            line_to_write += [gene_dict['gene_symbol'][i]]
            line_to_write += [gene_dict['locus'][i]]
            line_to_write += [gene_dict['TSS'][i]]
            line_to_write += [gene_dict['TES'][i]]
            line_to_write += [gene_dict['length'][i]]
            line_to_write += [str(nb_sample_above_10[i2])]
            line_to_write += [str(nb_sample_above_25[i2])]
            line_to_write += [str(nb_sample_above_50[i2])]
            line_to_write += [str(nb_sample_above_75[i2])]
            line_to_write += gene_dict['FPKM'][i]
            
            #print(line_to_write)
            fout.write('\t'.join(line_to_write)+'\n')
        
        fout.write('Total_FPKM\t\t\t\t\t\t\t'+str(total_above_1_fpkm)+'\t'+str(total_above_1_fpkm)+'\t'+str(total_above_1_fpkm)+'\t'+str(total_above_1_fpkm)+'\t'+'\t'.join([str(val) for val in sample_total])+'\n')
        fout.write('Main Isoform FPKM\t\t\t\t\t\t\t\t\t\t\t'+'\t'.join([str(val) for val in main1_fpkm])+'\n')
        if len(ids) > 1:
            fout.write('2nd Main Isoform FPKM\t\t\t\t\t\t\t\t\t\t\t'+'\t'.join([str(val) for val in main2_fpkm])+'\n')
        fout.write('Main Isoform/Total\t\t\t\t\t\t\t\t\t\t\t'+'\t'.join([str(val) for val in main1_pct])+'\n')
        if len(ids) > 1:
            fout.write('2nd Main Isoform / Main\t\t\t\t\t\t\t\t\t\t\t'+'\t'.join([str(val) for val in main2_over_main1])+'\n')
        fout.write('\n')

    fout.close()

if __name__ == '__main__':
	main()


