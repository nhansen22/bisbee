import pyensembl
import Bio.SeqIO
import Bio.Seq
import pandas as pd
import sys
import re
import math
import os
import io
from Bio import pairwise2
import utils as bb

print(sys.version_info)
print(pd._version)
print(bb.__file__)

events_file=sys.argv[1]
event_type=sys.argv[2]
aapad=int(sys.argv[3])
out_name=sys.argv[4]
ensemble_release=int(sys.argv[5])
ref_fasta=sys.argv[6]
species = "hg19"

#Functions from utils
def get_transcript_adj_exons(ensembl,gene_id,exon_coord):
 try:
  transcript_ids=ensembl.transcript_ids_of_gene_id(gene_id)
 except:
  print('Warning: ' + gene_id + ' not found')
  transcript_ids=[]
 transcript_list=[]
 for tid in transcript_ids:
   transcript=ensembl.transcript_by_id(tid)
   transcript.exon_intervals.sort()
   if exon_coord[0] in transcript.exon_intervals:
    idx=transcript.exon_intervals.index(exon_coord[0])
    if exon_coord == transcript.exon_intervals[idx:idx+len(exon_coord)]:
      transcript_list.append(transcript)
 return transcript_list

def has_coding_transcript(transcript_list):
 has_coding=False
 for transcript in transcript_list:
  if transcript.biotype=='protein_coding':
   has_coding=True
 return has_coding

def get_transcript_contain_exons(ensembl,gene_id,exon_coord):
 try:
  transcript_ids=ensembl.transcript_ids_of_gene_id(gene_id)
 except:
  print('Warning: ' + gene_id + ' not found')
  transcript_ids=[]
 transcript_list=[]
 for tid in transcript_ids:
   transcript=ensembl.transcript_by_id(tid)
   if set(exon_coord).issubset(set(transcript.exon_intervals)):
     transcript_list.append(transcript)
 return transcript_list

def find_overlap(rangeDF,coord):
 return (coord[0]<=rangeDF.loc[:,'end']) & (coord[1]>=rangeDF.loc[:,'start'])

def make_seq_from_coord(ref,contig,coordDF,strand):
  seq=''
  if not contig in ref:
   contig="chr"+str(contig)
   if not contig in ref:
    print(contig + "not found in ref")
    return seq
  for index,row in coordDF.iterrows():
   if strand=='+':
    seq=seq+str(ref[contig].seq[int(row.start)-1:int(row.end)])
   else:
    seq=str(ref[contig].seq[int(row.start)-1:int(row.end)].reverse_complement())+seq
  return seq

def find_seq_diff(seq1,seq2):
 align_list=pairwise2.align.globalms(seq1,seq2,2,-1,-10,0)
 if len(align_list)==0:
  seq1_diff_pos=(0,len(seq1)-1)
  seq2_diff_pos=(0,len(seq2)-1)
 elif seq1==seq2:
  seq1_diff_pos=(float('nan'),float('nan'))
  seq2_diff_pos=(float('nan'),float('nan'))
 else:
  align_data=pairwise2.format_alignment(*align_list[0]).split('\n')
  #print(align_data)
  letters= 'ACDEFGHIKLMNPQRSTVWY'
  seq_comp=pd.DataFrame({"seq1": list(align_data[0]), "seq2": list(align_data[2])})
  seq_comp=seq_comp.assign(seq1_pos= seq_comp.seq1.isin(list(letters)).cumsum())
  seq_comp=seq_comp.assign(seq2_pos= seq_comp.seq2.isin(list(letters)).cumsum())
  seq_comp=seq_comp.assign(match=seq_comp.seq1==seq_comp.seq2)
 #print(seq_comp.head())
 #print(seq_comp[seq_comp.match==False])
  first_mismatch=seq_comp.loc[seq_comp.match==False].index[0]
  if first_mismatch==0:
    seq1_diff_pos=(-1,max(seq_comp.seq1_pos[seq_comp.match==False]))
    seq2_diff_pos=(-1,max(seq_comp.seq2_pos[seq_comp.match==False]))
  else:
   seq1_diff_pos=(seq_comp.seq1_pos[first_mismatch-1],max(seq_comp.seq1_pos[seq_comp.match==False]))
   seq2_diff_pos=(seq_comp.seq2_pos[first_mismatch-1],max(seq_comp.seq2_pos[seq_comp.match==False]))

 return seq1_diff_pos, seq2_diff_pos

def is_coding_effect(transcript,effect_coord):
 coding_effect=dict()
 if transcript.biotype!='protein_coding' or not transcript.contains_stop_codon or not transcript.contains_start_codon:
  coding_effect['type']='NoncodingOrIncompleteTranscript'
  coding_effect['is_coding']=False
 else:
  coding_coord=pd.DataFrame(transcript.coding_sequence_position_ranges,columns=['start','end'])
  coding_effect['coord']=coding_coord.append(pd.Series({'start':min(transcript.stop_codon_positions),'end':max(transcript.stop_codon_positions)}),ignore_index=True).sort_values(by="start")
  coding_effect['overlap']=find_overlap(coding_effect['coord'],effect_coord)
  if effect_coord[1]<coding_coord.start.min() or effect_coord[0]>coding_coord.end.max():
   coding_effect['type']='UTR'
   coding_effect['is_coding']=False
  else:
   coding_effect['is_coding']=True
 return coding_effect

def make_prot_from_coord(transcript,coord,ref):
  trans=Bio.Seq.translate(make_seq_from_coord(ref,transcript.contig,coord,transcript.strand))
  stop_count=trans.count('*')
  if trans.startswith('M'):
   start_lost=False
  else:
   start_lost=True
  if stop_count==0:
   stop_lost=True
  else:
   stop_lost=False

  if stop_lost and not start_lost:
   if transcript.strand=='+':
    coord.iloc[-1,1]=transcript.exon_intervals[-1][1]
   else:
    coord.iloc[0,0]=transcript.exon_intervals[0][0]
    trans=Bio.Seq.translate(make_seq_from_coord(ref,transcript.contig,coord,transcript.strand))
    stop_count=trans.count('*')

  if start_lost or stop_count==0:
   prot_seq=''
  else:
   prot_seq=trans.split('*')[0]+'*'

  if start_lost:
   effect='StartLost'
  elif stop_count==0:
   effect='StopLost'
  elif stop_lost:
   effect='PostStop'
  elif stop_count==1 and trans.endswith('*'):
   effect='InFrame'
  else:
   effect='PrematureStop'
  return prot_seq, effect

def get_event_coords(event_info,event_type):
 event_coords=pd.DataFrame(columns=("isoform","start","end"))
 if event_type=="IR":
  event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon1_end"],"end":event_info["exon2_start"]},ignore_index=True)
 elif event_type=="ES":
  event_coords=event_coords.append({"isoform":"iso1","start":event_info["exon_pre_end"],"end":event_info["exon_aft_start"]},ignore_index=True)
  event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon_pre_end"],"end":event_info["exon_start"]},ignore_index=True)
  event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon_end"],"end":event_info["exon_aft_start"]},ignore_index=True)
 elif event_type=="MUT":
  event_coords=event_coords.append({"isoform":"iso1","start":event_info["exon_pre_end"],"end":event_info["exon1_start"]},ignore_index=True)
  event_coords=event_coords.append({"isoform":"iso1","start":event_info["exon1_end"],"end":event_info["exon_aft_start"]},ignore_index=True)
  event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon_pre_end"],"end":event_info["exon2_start"]},ignore_index=True)
  event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon2_end"],"end":event_info["exon_aft_start"]},ignore_index=True)
 elif (event_type=="A3" and event_info["strand"]=="+") or (event_type=="A5" and event_info["strand"]=="-"):
  event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon_const_end"],"end":event_info["exon_alt1_start"]},ignore_index=True)
  event_coords=event_coords.append({"isoform":"iso1","start":event_info["exon_const_end"],"end":event_info["exon_alt2_start"]},ignore_index=True)
 elif (event_type=="A3" and event_info["strand"]=="-") or (event_type=="A5" and event_info["strand"]=="+"):
  event_coords=event_coords.append({"isoform":"iso2","start":event_info["exon_alt1_end"],"end":event_info["exon_const_start"]},ignore_index=True)
  event_coords=event_coords.append({"isoform":"iso1","start":event_info["exon_alt2_end"],"end":event_info["exon_const_start"]},ignore_index=True)
 return event_coords

def jid_to_coords(event_jid):
 event_coords=pd.DataFrame(columns=("isoform","start","end"))
 iso1_coords=event_jid.split("g.")[1].split(">")[0].split('_')
 for coord in iso1_coords:
  if not coord=="NONE":
   event_coords=event_coords.append({"isoform":"iso1","start":int(coord.split('j')[0]),"end":int(coord.split('j')[1])},ignore_index=True)
 iso2_coords=event_jid.split("g.")[1].split(">")[1].split("[")[0].split('_')
 for coord in iso2_coords:
  if not coord=="NONE":
   event_coords=event_coords.append({"isoform":"iso2","start":int(coord.split('j')[0]),"end":int(coord.split('j')[1])},ignore_index=True)
 return event_coords


def find_matching_transcripts(ensembl,gene_id,event_coords):
 try:
  transcript_ids=ensembl.transcript_ids_of_gene_id(gene_id[0:15])
 except:
  print('Warning: ' + gene_id[0:15] + ' not found')
  transcript_ids=[]
 transcript_table=pd.DataFrame(columns=["coding","matching_isoform"],index=transcript_ids)
 for tid in transcript_ids:
  transcript=ensembl.transcript_by_id(tid)
  if transcript.biotype!='protein_coding' or not transcript.contains_stop_codon or not transcript.contains_start_codon:
   transcript_table.loc[tid,"coding"]=False
  else:
   transcript_table.loc[tid,"coding"]=True
  transcript_table.loc[tid,"matching_isoform"]=get_matching_isoform(transcript,event_coords)
 return transcript_table

def get_matching_isoform(transcript,event_coords):
 exons=pd.DataFrame(transcript.exon_intervals,columns=["start","end"]).sort_values(by="start")
 event_region=[event_coords.start.min(),event_coords.end.max()]
 exons=exons[find_overlap(exons,event_region)]
 junctions=pd.DataFrame(columns=["start","end"])
 junctions["end"]=exons.start[1:].values
 junctions["start"]=exons.end[0:-1].values
 if junctions.equals(event_coords.loc[event_coords.isoform=="iso2",["start","end"]].reset_index(drop=True).astype(int)):
  matching_isoform="iso2"
 elif sum(event_coords.isoform=="iso1")==0:
  if len(exons.start)>0:
   if event_coords.start[0]>exons.start.iloc[0] and event_coords.end[0]<exons.end.iloc[0]:
    matching_isoform="iso1"
   else:
    matching_isoform="none"
  else:
    matching_isoform="none"
 elif junctions.equals(event_coords.loc[event_coords.isoform=="iso1",["start","end"]].reset_index(drop=True).astype(int)):
   matching_isoform="iso1"
 else:
   matching_isoform="none"
 return matching_isoform

def get_new_coord(matching_isoform,event_coords,old_coord):
 novel_junc=event_coords.loc[event_coords.isoform!=matching_isoform]
 overlap=find_overlap(old_coord,[event_coords.start.min(),event_coords.end.max()])
 if novel_junc.size>0 and overlap.any():
  new_coord=pd.DataFrame(columns=["start","end"])
  new_coord=new_coord.append({'start':old_coord[overlap].start.iloc[0],'end':novel_junc.start.iloc[0]},ignore_index=True)
  new_coord=new_coord.append(pd.DataFrame({'start':novel_junc.end.iloc[0:-1].values,'end':novel_junc.start.iloc[1:].values}),ignore_index=True)
  new_coord=new_coord.append({'start':novel_junc.end.iloc[-1],'end':old_coord[overlap].end.iloc[-1]},ignore_index=True)
  new_coord=new_coord.append(old_coord.loc[overlap==False],ignore_index=True).sort_values(by="start")
 elif overlap.any():
  new_coord=pd.DataFrame(columns=["start","end"])
  new_coord=new_coord.append({'start':old_coord[overlap].start.iloc[0],'end':old_coord[overlap].end.iloc[-1]},ignore_index=True)
  new_coord=new_coord.append(old_coord.loc[overlap==False],ignore_index=True).sort_values(by="start")
 else:
  new_coord=old_coord
 return new_coord

if species=="hg19":
    ensembl=pyensembl.genome.Genome("GRCh37", "ensembl74", annotation_version=74, gtf_path_or_url="ftp://ftp.ensembl.org/pub/release-74/gtf/homo_sapiens/Homo_sapiens.GRCh37.74.gtf.gz", transcript_fasta_paths_or_urls="/home/nhansen/Packages/pyensembl/Homo_sapiens.GRCh37.74.cdna.all.fa.gz", protein_fasta_paths_or_urls="/home/nhansen/Packages/pyensembl/Homo_sapiens.GRCh37.74.pep.all.fa.gz", decompress_on_download=False, copy_local_files_to_cache=False, cache_directory_path=None)
elif species=="hg38":
    ensembl=pyensembl.genome.Genome("GRCh38", "ensembl98", annotation_version=98, gtf_path_or_url="/home/nhansen/.cache/pyensembl/GRCh38/ensembl98/Homo_sapiens.GRCh38.98.gtf.gz", transcript_fasta_paths_or_urls="/home/nhansen/.cache/pyensembl/GRCh38/ensembl98/Homo_sapiens.GRCh38.cdna.all.fa.gz", protein_fasta_paths_or_urls="/home/nhansen/.cache/pyensembl/GRCh38/ensembl98/Homo_sapiens.GRCh38.pep.all.fa.gz", decompress_on_download=False, copy_local_files_to_cache=False, cache_directory_path=None)
else:
    print("ERROR: invalid species selection, choose mm10, hg19, hg38!!!")
    
#ensembl=pyensembl.EnsemblRelease(ensemble_release)
ref=Bio.SeqIO.to_dict(Bio.SeqIO.parse(ref_fasta,'fasta'))

###read coordinates table
if event_type=="IR":
 col_types={"contig":str,"event_id":str,"gene":str,"exon1_end":int,"exon2_start":int,"event_jid":str}
elif event_type=="ES":
 col_types={"contig":str,"event_id":str,"gene":str,"exon_pre_end":int,"exon_start":int,"exon_end":int,"exon_aft_start":int,"event_jid":str}
elif event_type=="MUT":
 col_types={"contig":str,"event_id":str,"gene":str,"exon_pre_end":int,"exon1_start":int,"exon1_end":int,"exon2_start":int,"exon2_end":int,"exon_aft_start":int,"event_jid":str}
elif event_type=="A3" or event_type=="A5":
 col_types={"contig":str,"event_id":str,"gene":str,"strand":str,"exon_const_start":int,"exon_const_end":int,"exon_alt1_start":int,"exon_alt1_end":int,"exon_alt2_start":int,"exon_alt2_end":int,"event_jid":str}
elif event_type=="ALL":
 col_types={"contig":str,"event_id":str,"gene":str,"strand":str,"event_jid":str}
else:
 print("ERROR: invalid event type " + event_type + ".\nValid event types are IR, ES, MUT, A3 or A5")
 sys.exit(1)

events_table=pd.read_csv(events_file,usecols=col_types.keys(),dtype=col_types)
events_table=events_table.assign(coding_transcript_effect=None,top_effect_transcript=None,effect_type=None,iso1_pc=None,iso2_pc=None,other_pc=None,iso1_nc=None,iso2_nc=None,other_nc=None)
topEffect_list=[]
wt_file='.'.join([out_name,event_type,"refSeq.fasta"])
if os.path.exists(wt_file):
 os.remove(wt_file)
wt_fasta=open(wt_file,"a+")
novel_file='.'.join([out_name,event_type,"altSeq.fasta"])
if os.path.exists(novel_file):
 os.remove(novel_file)
novel_fasta=open(novel_file,"a+")
peptides_file='.'.join([out_name,event_type,"peptides.csv"])
if os.path.exists(peptides_file):
 os.remove(peptides_file)
peptides_table=open(peptides_file,"a+")
effectCol=["altPept","event_id","effectId","gene","orf_effect","aa_effect","topEffect","refPept","refIsoform","sourceName","insSeqLen","refSeqLen","delSeqLen","refSeqPos","altSeqPos"]
effectDF=pd.DataFrame(columns=effectCol)
effectDF.to_csv(peptides_table,header=True,index=False)
 
 
## for each event
###find matching transcripts
###determine effect_type
for index,row in events_table.iterrows():
 effectDF=pd.DataFrame(columns=effectCol)
 wt_list=[]
 novel_list=[]
 event_coords=jid_to_coords(events_table.loc[index,"event_jid"])
 #event_coords=bb.get_event_coords(row,event_type)
 iso1_str=events_table.loc[index,"event_jid"].split("g.")[1].split(">")[0]
 iso2_str=events_table.loc[index,"event_jid"].split("g.")[1].split(">")[1].split("[")[0]
 transcript_table=find_matching_transcripts(ensembl,row["gene"],event_coords)
 events_table.loc[index,"iso1_pc"]='|'.join(transcript_table.index.values[(transcript_table.coding==True) & (transcript_table.matching_isoform=="iso1")])
 events_table.loc[index,"iso2_pc"]='|'.join(transcript_table.index.values[(transcript_table.coding==True) & (transcript_table.matching_isoform=="iso2")])
 events_table.loc[index,"other_pc"]='|'.join(transcript_table.index.values[(transcript_table.coding==True) & (transcript_table.matching_isoform=="none")])
 events_table.loc[index,"iso1_nc"]='|'.join(transcript_table.index.values[(transcript_table.coding==False) & (transcript_table.matching_isoform=="iso1")])
 events_table.loc[index,"iso2_nc"]='|'.join(transcript_table.index.values[(transcript_table.coding==False) & (transcript_table.matching_isoform=="iso2")])
 events_table.loc[index,"other_nc"]='|'.join(transcript_table.index.values[(transcript_table.coding==False) & (transcript_table.matching_isoform=="none")])
 if events_table.loc[index,"iso1_pc"]!='' and events_table.loc[index,"iso2_pc"]!='':
  events_table.loc[index,"coding_transcript_effect"]="IsoformSwitch"
 elif events_table.loc[index,"iso1_pc"]!='':
  events_table.loc[index,"coding_transcript_effect"]="NovelIso2"
 elif events_table.loc[index,"iso2_pc"]!='':
  events_table.loc[index,"coding_transcript_effect"]="NovelIso1"
 elif events_table.loc[index,"other_pc"]!='':
  events_table.loc[index,"coding_transcript_effect"]="Unknown"
 else:
  events_table.loc[index,"coding_transcript_effect"]="NonCoding"
 #print(events_table.loc[index,])
 top_coding_effect=''
 top_effectId=''
 top_effect_type=''
 top_effect_score=-5
 max_wt_seq=-1
 for tid,transcript_info in transcript_table.loc[(transcript_table.coding==True) & ((transcript_table.matching_isoform=="iso1") | (transcript_table.matching_isoform=="iso2"))].iterrows():
  #print(tid)
  curr_effect_score=-5
  transcript=ensembl.transcript_by_id(tid)
  wt_seq=Bio.Seq.translate(transcript.coding_sequence)
  coding_coord=pd.DataFrame(transcript.coding_sequence_position_ranges,columns=['start','end'])
  coding_coord=coding_coord.append(pd.Series({'start':min(transcript.stop_codon_positions),'end':max(transcript.stop_codon_positions)}),ignore_index=True).sort_values(by="start")
  new_coord=get_new_coord(transcript_info["matching_isoform"],event_coords,coding_coord)
   #print(novel_junc)
   #print(coding_coord)
   #print(overlap)
   #print(new_coord)
  (mut_seq,coding_effect)=make_prot_from_coord(transcript,new_coord,ref)
   #print(mut_seq)
  (wt_diff_pos,mut_diff_pos)=find_seq_diff(wt_seq,mut_seq)
   #print(wt_diff_pos)
   #print(mut_diff_pos)
  desc=' '.join(["GN=" + transcript.gene.gene_name,"ID=" + events_table.loc[index,"event_id"] + "_" + transcript_info["matching_isoform"],"JID=" + events_table.loc[index,"event_jid"],"PC=" + str(wt_diff_pos[0]) + "-" + str(wt_diff_pos[1])])
  wt_list.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id,description=desc,seq=Bio.Seq.Seq(wt_seq)))
  if transcript_info["matching_isoform"]=="iso1":
   desc=' '.join(["GN=" + transcript.gene.gene_name,"ID=" + events_table.loc[index,"event_id"] + "_iso2","JID=" + events_table.loc[index,"event_jid"],"PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])
   novel_list.append(Bio.SeqRecord.SeqRecord(id= "en|" +transcript.id + ":g." + iso1_str + ">" + iso2_str,description=desc,seq=Bio.Seq.Seq(mut_seq)))
  else:
   desc=' '.join(["GN=" + transcript.gene.gene_name,"ID=" + events_table.loc[index,"event_id"] + "_iso1","JID=" + events_table.loc[index,"event_jid"],"PC=" + str(mut_diff_pos[0]) + "-" + str(mut_diff_pos[1])])
   novel_list.append(Bio.SeqRecord.SeqRecord(id="en|" + transcript.id + ":g." + iso2_str + ">" + iso1_str,description=desc,seq=Bio.Seq.Seq(mut_seq)))
#"mutPept","event_id","event_jid","gene_name","transcript_id","orf_effect","topEffect","wtPept","sourceName","novelSeqLen","wtSeqLen","delSeqLen"]
  if not math.isnan(mut_diff_pos[0]):
   mut_pept=mut_seq[max(mut_diff_pos[0]-aapad,0):mut_diff_pos[1]+aapad]
   wt_pept=wt_seq[max(mut_diff_pos[0]-aapad,0):wt_diff_pos[1]+aapad]
   novel_seq_len=mut_diff_pos[1]-mut_diff_pos[0]
   del_seq_len=wt_diff_pos[1]-wt_diff_pos[0]
  else:
   mut_pept=''
   wt_pept=''
   novel_seq_len=float('nan')
   del_seq_len=float('nan')
  if (coding_effect=="InFrame") and math.isnan(mut_diff_pos[0]):
   aa_effect="Silent"
   curr_effect_score=-4
  elif del_seq_len==len(wt_seq)-1:
   aa_effect="ProteinLoss"
   curr_effect_score=-2
  elif (novel_seq_len>0) and (del_seq_len>0):
   aa_effect="Substitution"
   curr_effect_score=novel_seq_len
  elif (novel_seq_len>0):
   aa_effect="Insertion"
   curr_effect_score=novel_seq_len
  elif (del_seq_len>0) and (wt_diff_pos[1]==len(wt_seq)-1):
   aa_effect="Truncation"
   curr_effect_score=-3
  elif (del_seq_len>0):
   aa_effect="Deletion"
   curr_effect_score=-1
  else:
   aa_effect="Unknown"
  effectId=tid + "_" +  events_table.loc[index,"event_jid"]
  if (curr_effect_score>top_effect_score) or ((curr_effect_score==top_effect_score) and (max_wt_seq>len(wt_seq))):
   top_coding_effect=aa_effect + "_" + coding_effect
   top_effectId=effectId
   top_effect_score=curr_effect_score
   max_wt_seq=len(wt_seq)
   curr_effect=pd.Series({"altPept": mut_pept, "event_id": events_table.loc[index,"event_id"], "effectId": effectId, "gene": transcript.gene.gene_name,"orf_effect": coding_effect,"aa_effect": aa_effect,"refPept": wt_pept,"refIsoform": transcript_info["matching_isoform"],"insSeqLen": novel_seq_len,"refSeqLen": len(wt_seq)-1,"delSeqLen":del_seq_len,"refSeqPos":str(wt_diff_pos[0]+1) + '-' + str(wt_diff_pos[1]+1),"altSeqPos":str(mut_diff_pos[0]+1) + '-' + str(mut_diff_pos[1]+1)})
   #print(curr_effect)
   effectDF=effectDF.append(curr_effect,ignore_index=True)

 out_handle=io.StringIO()
 Bio.SeqIO.write(wt_list,out_handle,"fasta")
 wt_fasta.write(out_handle.getvalue())
 out_handle=io.StringIO()
 Bio.SeqIO.write(novel_list,out_handle,"fasta")
 novel_fasta.write(out_handle.getvalue())
 effectDF.topEffect=effectDF.effectId==top_effectId
 effectDF.to_csv(peptides_table,header=False,index=False)
 events_table.loc[index,"effect_type"]=top_coding_effect
 events_table.loc[index,"top_effect_transcript"]=top_effectId.split('_')[0]


events_table.to_csv('.'.join([out_name,event_type,'effects.csv']))
wt_fasta.close()
novel_fasta.close()
peptides_table.close()
#Bio.SeqIO.write(wt_list,'.'.join([out_name,event_type,"wtSeq.fasta"]),"fasta")
#Bio.SeqIO.write(novel_list,'.'.join([out_name,event_type,"novelSeq.fasta"]),"fasta")
#### for each isoform switch
#####add seq to fastas

#### for each iso1/iso2 seq
#####make new coord
#####make new seq
#####make effect id
#####compare seq
#####classify effect
#####add seq to fastas
#####get mut/wt peptides
#####score effect

### find top effect

