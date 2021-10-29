from Bio import SeqIO
from time import sleep
from Bio import AlignIO, Align
from Bio.Align.AlignInfo import PSSM, SummaryInfo
from io import StringIO
from Bio import Entrez
from tqdm import tqdm

# Unix systems ONLY
from sh import mafft

import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import logging 
import time

logger = logging.getLogger()

class FeatureFinder(object):
    """Feature method to read GenBank File Formats into a dataframe"""
    def __init__(self):
        self.rna_df = None

    def get_features(self, gbff_file_path, drop_duplicates=True):
        """Primary function to feed the gbff file path to."""
        rna_df = self.get_gbff_features(gbff_file_path, drop_duplicates=drop_duplicates)
        rna_df = self.get_df_lengths(rna_df)
        self.rna_df = rna_df
        return rna_df

    def extract_cds(self, gb_feature, gb_recs):
        # get all the annotated qualifiers within the gb records
        qualifiers = list(gb_feature.qualifiers.keys())

        if 'gene' in qualifiers:
            cds = ','.join(gb_feature.qualifiers['gene'])
        elif 'product' in qualifiers:
            cds = gb_feature.qualifiers['product']
            cds = ','.join([cd.replace(' ', '_') for cd in cds])
        else:
            cds = gb_feature.type

        if 'locus_tag' in qualifiers:
            locus = ','.join(gb_feature.qualifiers['locus_tag'])
        else:
            locus = None

        # location of the coding sequence within the genome
        cds_loc = str(gb_feature.location)
        if 'join' in cds_loc:
            return None
        cds_start = gb_feature.location.start.position
        cds_end = gb_feature.location.end.position
        # parse the record to get the coding piece
        cds_seq = gb_recs[cds_start:cds_end]
        # check the orientation of the cds

        if gb_feature.strand == -1:
            orientation = 'reverse'
            cds_seq = cds_seq.reverse_complement().seq.__str__()
        else:
            orientation = 'forward'
            cds_seq = cds_seq.seq.__str__()

        cds_data = [cds, locus, cds_loc, cds_start, cds_end, orientation, cds_seq]
        return cds_data
    
    def get_df_lengths(self, rna_df):
        rna_df['length'] = rna_df.End - rna_df.Start
        rna_df['aa_length'] = rna_df.length / 3
        rna_df2 = rna_df[(rna_df.Gene.str.contains('ypothetical') == False) & (rna_df.Gene != 'gene')]
        rna_df2 = rna_df2.sort_values(by='aa_length')

        return rna_df2

    def get_gbff_features(self, gbff_handle, drop_duplicates=True):

        # use the next iterator to load all the records within the genbank file
        gb_master = SeqIO.parse(gbff_handle, 'genbank')
        # define a function to extract a CDS.

        # assemble all the genes into a list
        genes = []
        types = []
        # loop over all features in the loaded gb records
        for gb_recs in gb_master:
            for feature in gb_recs.features:
                if feature.type in ['CDS', 'ncRNA', 'rRNA', 'miRNA', 'mRNA', 'gene']:
                    # this is a coding sequence
                    feat = self.extract_cds(feature, gb_recs)
                    # if this is not an empty list, add
                    if feat:
                        genes.append(feat)
                        types.append(feature.type)

        # put into a df which enumerates all the coding transcripts
        df_transcripts = pd.DataFrame(genes, columns=['Gene', 'Locus', 'Location',
                                                      'Start', 'End', 'Orientation', 'Sequence'])
        df_transcripts['Feature'] = types
        if drop_duplicates is True:
            df_transcripts.drop_duplicates(subset=['Gene', 'Locus'], keep='last', inplace=True)
        return df_transcripts

class EntrezWrapper(object):
    """
    Downloads files from the NCBI nucleotide database.
    Input:
    email: str, entrez requires an email for every query
    ncbi_api_key: str, strongly reccomended to include this for large queries
    """
    def __init__(self, email, ncbi_api_key=None):
        self.entrez = Entrez
        self.entrez.email = email
        if ncbi_api_key is not None:
            self.entrez.api_key = ncbi_api_key
        self.entrez.max_tries = 3
        self.entrez.sleep_between_tries = 3

    def get_search_dict(self, database, mquery, sort='significance'):
        handle = self.entrez.esearch(database, term=mquery, usehistory="y", sort=sort)
        search_results = Entrez.read(handle)
        s_count = int(search_results['Count'])
        message = f'Found {s_count} entries associated with {mquery} in {database} database'
        print(message)
        return search_results, s_count

    def get_nucleotide_fasta(self, mquery, outfile, db='nucleotide', write=False, retmax=300, limit=300):
            """"
            Input:
            query: str, ncbi nucleotide search term, should be pathogen
            outfile: str, path to output fasta file
            write: bool, write output to file path otherwise return as string.
            """

            search_results, s_count = self.get_search_dict(db, mquery)
            time.sleep(1)

            if s_count == 0:
                message = f'Query invalid: {query}'
                logger.error(message)
                raise ValueError(message)

            master_fetch = ''
            for start in range(0, limit, retmax):
                print(start, 'quer_lim', retmax, limit)
                if start >= limit:
                    break
                fetch_handle = self.entrez.efetch(db=db, rettype="fasta", retmode="text", retstart=start, retmax=retmax,
                                                  webenv=search_results['WebEnv'], query_key=search_results['QueryKey'])
                time.sleep(1)
                master_fetch += fetch_handle.read()
                master_fetch += '\n'
            if write is True:
                with open(outfile, 'w') as fasta_file:
                    fasta_file.write(master_fetch)
            else:
                return master_fetch
            
class PullQCGenes(object):
    def __init__(self, email='None@gmail.com', ncbi_api_key=None):
        self.df_list = []
        self.rec_count = 0
        self.ee = EntrezWrapper(email, ncbi_api_key=ncbi_api_key)
        self.completed_genes = {}

    @staticmethod
    def remove_seq_duplicates(fastas, seq_len, format='fasta'):
        fasta_array = [rec for rec in SeqIO.parse(StringIO(fastas), 'fasta') if len(rec) < seq_len * 1.5]
        fasta_df = pd.DataFrame([([rec.description, rec.seq.__str__()]) for rec in fasta_array], columns=['description', 'seq'])
        fasta_df.drop_duplicates(subset='seq', inplace=True)
        if format == 'fasta':
            fastas2 = "\n>".join(fasta_df['description'] + '\n' + fasta_df['seq'])
            return fastas2
        else:
            return fasta_df
        
    def plot_pssm(self, df_pssm, protein_name):
        plt.figure(figsize=(30,4),dpi=300)
        sns.heatmap(df_pssm)
        plt.ylabel('AA at position')
        plt.xlabel('Position in alignment')
        plt.title(protein_name)
        plt.show()
        
    def pull_and_qc_genes(self, rna_df2, organism, retmax=300, limit=300):
        ee = self.ee
        df_list = self.df_list
    
        for row in rna_df2.itertuples():
            protein_name = row[1]
            seq_len = row[-1]
            if protein_name in self.completed_genes:
                continue
            else:
                self.completed_genes[protein_name] = None
            try:
                fastas = ee.get_nucleotide_fasta(f'{protein_name} AND {organism}[Organism]', f'{protein_name}_seq.fasta', db='protein', 
                                                write=False, retmax=retmax, limit=limit)
            except Exception as e:
                print("Exception", e)
                continue

            # cut down on length
            fasta_array = [rec for rec in SeqIO.parse(StringIO(fastas), 'fasta') if len(rec) < seq_len * 1.5 and len(rec) > seq_len *0.5]
            fastas = "\n".join([">" + "\n".join([rec.description, rec.seq.__str__()]) for rec in fasta_array])

            fastas = self.remove_seq_duplicates(fastas, seq_len)

            filepath = f'{protein_name}.fasta'
            open(filepath, 'w').write(fastas)
            buf = StringIO(fastas)
            bufout = StringIO()

            # wrap mafft alignment
            mafft(filepath, _out=bufout)
            bufout.seek(0)

            aln = bufout.read()
            aln_df = self.remove_seq_duplicates(aln, seq_len, format='df')
            aln_df['gene'] = protein_name
            bufout.seek(0)
            try:
                alnio = AlignIO.read(bufout,'fasta')
            except Exception as e:
                print(e)
                continue
            pssm = (SummaryInfo(alnio).pos_specific_score_matrix().pssm)
            pssm2 = {k+str(ind):v for ind, (k,v) in enumerate(pssm)}
            df_pssm = pd.DataFrame(pssm2)

            self.plot_pssm(df_pssm, protein_name)
            row_len, col_len = df_pssm.shape

            if col_len < seq_len * 2 and row_len > 4:
                self.rec_count += row_len
                df_list.append(aln_df)
            else:
                print(f'{protein_name} retrieved information insufficient')

            print("Total rec count: ", self.rec_count)
            
        df = pd.concat(df_list)
        
        return df
 
