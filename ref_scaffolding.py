#!/usr/bin/env python 
import os,sys
import argparse 

class HelpFormatter(argparse.RawDescriptionHelpFormatter,argparse.ArgumentDefaultsHelpFormatter):
    pass

# MUMmer version >= 4.0.0 include nucmer show-coords *  
_MUMMER = "/export/personal/fanjp/software/MUMmer/current/bin"
_AUTHER = "Minghui Meng"
_MAIL   = "hyacinlee@163.com"

def main():
    
    args  = get_args()
    ref   = os.path.abspath(args.r)
    query = os.path.abspath(args.a)
    outdir= os.path.abspath(args.o)

    if not os.path.exists(outdir):
        os.mkdir(outdir)
    
    if 1 in args.s:
        align(ref,query,outdir+"/1.align/",args.t)

    if 2 in args.s:
        creat_new_genome(query,ref,outdir)

    if 3 in args.s:
        align(ref,outdir+"/2.scaffolding/assemblyed.chroms.fasta",outdir+"/3.align/",args.t)

    if 4 in args.s:
        draw_dotplot(outdir+"/3.align/",args.f,outdir+"/2.scaffolding/assemblyed.chroms.fasta")



def draw_dotplot(outdir,filter,fasta):
    print "# Draw dotplot between two geneome ... ..."
    run_or_die("%s/delta-filter -l %s -i 90 %s/align.delta > %s/align.f%s.delta"    % (_MUMMER, filter , outdir ,outdir ,filter) )
    run_or_die("cd %s && %s/mummerplot -f --png -p align.f%s %s/align.f%s.delta" % (outdir, _MUMMER, filter , outdir ,filter) )
    for chromsome in read_fasta(fasta,"list"):
        run_or_die("cd %s && %s/mummerplot -f --png -r %s -q %s -p align.f%s.%s  %s/align.f%s.delta" % (outdir, _MUMMER, chromsome, chromsome, filter, chromsome , outdir ,filter) )



def creat_new_genome(query,ref,outdir):

    if not os.path.exists(outdir+"/2.scaffolding/"):
        os.mkdir(outdir+"/2.scaffolding/")

    print "# Reading align information to find contig's best match chromsome ... ... "
    (info,best) = read_align_result(outdir+"/1.align/align.coords")

    print "# Reading align information to find contig's best match location and direction  ... ..."
    order = re_order(info,best) 

    print "# Reading chromsome order in reference geneome  ... ..."
    chr_order = read_fasta(ref,"list")

    (fasta,genome_size) = read_fasta(query,"dict")

    print "# Creat new genome ... ...."
    out1 = open(outdir+"/2.scaffolding/assemblyed.chroms.fasta","w")
    out2 = open(outdir+"/2.scaffolding/assemblyed.chroms.agp","w")
    out3 = open(outdir+"/2.scaffolding/assemblyed.chroms.stat","w")

    out3.write("#Chr\tNum\tlength\tlength(withN)\n")
    (tl_withoutN,tl_withN,tl_contig) = (0,0,0)

    for chromsome in chr_order:
        
        if chromsome not in order:
            continue

        (seq,chr_len,end) = ("",0,0)
        for i,info in enumerate(order[chromsome]):
            contig_seq = "".join(fasta[info[0]])
            start      = len(seq)+1

            if info[1] == "+":
                seq += contig_seq
            else:
                seq += reverse_fasta(contig_seq)

            end = len(seq) 
            if i < len(order[chromsome]) - 1 :
                seq += "N"*100
            
            chr_len     += len(contig_seq)
            tl_withoutN += len(contig_seq)
        
            out2.write("%s\t%s\t%s\t%s\t1\t%s\t%s\n" % (chromsome,start,end,info[0],len(contig_seq),info[1]) )
        
        tl_withN += len(seq)
        tl_contig   += len(order[chromsome])
        out1.write(">%s\n%s\n" %(chromsome,seq))
        out3.write("%s\t%s\t%s\t%s\n" % (chromsome,len(order[chromsome]),chr_len,end))

    out3.write("Total\t%s\t%s\t%s\n" % (tl_contig,tl_withoutN,tl_withN))
    out3.write("Load Rate\t%s\n" %(tl_withoutN/float(genome_size)))

    out1.close()
    out2.close()
    out3.close()



def reverse_fasta(sequence):

    transline=sequence[::-1].upper().replace('A','t').replace('T','a').replace('G','c').replace('C','g').upper()
    return transline


def re_order(info,best):

    postion = {}
    order   = {}

    for contig in info.keys():
        (locate,number) = (0,0)
        stand = {"+":0,"-":0}
        
        for align in info[contig]:
            if align[0] != best[contig]:    # not mapped to the best chromsome
                continue
            locate          += align[1]
            stand[align[3]] += 1

        locate = locate/(stand["+"]+stand["-"])
        st     = sorted(stand.items(),key=lambda x:x[1],reverse=True)[0][0]

        if best[contig] not in postion:
            postion[best[contig]] = [(contig,st,locate)]
        else:
            postion[best[contig]].append((contig,st,locate))

    for chromsome in postion.keys():
        sort_data = sorted(postion[chromsome],key=lambda x:(int(x[2])))
        order[chromsome] = [ (x[0],x[1]) for x in sort_data]

    return order



def read_align_result(coords_file):

    leng = {}
    info = {}
    best = {}

    for line in open(coords_file,"r"):

        lst = line.strip().replace("|","").split()
        if len(lst) != 13:
            continue

        length = int(lst[3]) -int(lst[2])
        stand  = "+"
        if length < 0:
            length = -length
            stand = "-"

        if lst[12] not in info:
            info[lst[12]] = [(lst[11],int(lst[0]),length,stand)]  # chr  , start , length , stand 
            leng.update( {lst[12]:{lst[11]:length}})
        else:
            info[lst[12]].append((lst[11],int(lst[0]),length,stand)) 

            if lst[11] not in leng[lst[12]]:
                leng[lst[12]].update({lst[11]:length})
            else:
                leng[lst[12]][lst[11]] += length

    # find best align chromsome for each contig         
    for key in leng.keys():   
        best[key] = sorted(leng[key].items(),key=lambda x:x[1],reverse=True)[0][0]
        #print key,best[key],leng[key]

    return info,best



def align(ref,query,outdir,t):

    os.system("mkdir -p "+outdir)
    cmd = "cd %s && %s/nucmer -g 1000 --prefix align -c 1000 -l 500 -t %s %s %s " %(outdir,_MUMMER,t,ref,query)
    run_or_die(cmd)
    cmd = "%s/show-coords -clr  -L 5000  %s/align.delta > %s/align.coords " %(_MUMMER,outdir,outdir)
    run_or_die(cmd)



def read_fasta(fasta_file,tp):

    print "# Reading fasta dict from %s" %(fasta_file)          # append to list
    chr_lst  = []
    fasta    = {}
    fasta_id = ''
    size     = 0
    for line in open(fasta_file):
        if line.startswith(">"):
            fasta_id = line.strip().replace(">","").split()[0]
            chr_lst.append(fasta_id)
            fasta[fasta_id] = []
        else:
            fasta[fasta_id].append(line.strip())
            size += len(line.strip())

    if tp == "list":
        return chr_lst
    elif tp == "dict":
        return (fasta,size)
    else:
        return (fasta,chr_lst)



def run_or_die(cmd):

    print "# Ready to Run: "+cmd

    flag = os.system(cmd)

    if flag != 0:
        print "# Error cmd: %s\n" %(cmd)
        exit(1)
    else:
        print "# done! "


def get_args():

    parser = argparse.ArgumentParser(
    formatter_class = HelpFormatter,
    description = '''
RunCmd :
    Funciton: Scaffoldding your assemblyed contig/scaffold genome base on reference align 
'''
    )
    parser.add_argument('-a',metavar='',help='assemblyed contig/scaffold genome')
    parser.add_argument('-r',metavar='',help='reference chromsomes genome')
    parser.add_argument('-o',metavar='',help='outdir',default="./output")
    parser.add_argument('-t',metavar='',help='thread of genome align by nucmer',default=20)
    parser.add_argument('-f',metavar='',help='min size to draw dotplot ',default=10000)
    parser.add_argument('-s',metavar='',help='step choice [1 -> align, 2 -> creat new genome, 3 -> align, 4 -> dotplot ]',nargs="+",default=[1,2,3,4])
    args = parser.parse_args()
    
    if not args.a or not args.r:
        parser.print_help()
        exit(1)

    return args


if __name__ == '__main__':
    main()