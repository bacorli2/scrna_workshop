

sra_id="SRR23691690"
transcriptome_ref_path="/home/bacorliss/research/zebra_scrna/ref/lawson/Danio-rerio-genome-lawson"

# This script downloads the specified SRA file, splits it into separate reads, 
# and then performs alignment with a zebrafish reference genome. All 
# processing is done in the current working directory, which is used as a temp 
# data folder. The output count matrix is exported to the parent directory, 
# and then all the downloaded and intermediate data is deleted.
# Progress is tracked with empty status files that are written to [sra_id]_status
# so (in theory) the script will continue where it left off if the processes is
# terminated.

# MIT License
#
# Copyright (c) [2023] [Bruce Allen Corliss, North Carolina State University]
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
mkdir ${sra_id}_status

# Fetch files from SRA
if ! test -f ./${sra_id}_status/1_prefectch_complete; then
  # Delete incomplete files if they exist
  if test -f ./${sra_id}; then
    rm -r ./${sra_id}
  fi

  # download SRA data
  prefetch -v ${sra_id} --max-size 500000000
  
  touch ./${sra_id}_status/1_prefectch_complete
fi


# Separate read files
if ! test -f ./${sra_id}_status/2_fastq-dump_complete; then
  # Delete incomplete files if they exist
  if test -f ${sra_id}_1.fastq; then
    rm -f ./${sra_id}_{1,4}.fastq
  fi
  
  fasterq-dump ./${sra_id} --split-files --include-technical
  
  touch ./${sra_id}_status/2_fastq-dump_complete
fi
#> ${sra_id}_1.fastq, ${sra_id}_2.fastq, 
#> ${sra_id}_3.fastq, ${sra_id}_4.fastq


# Rename read files for cell ranger
if ! test -f ./${sra_id}_status/3_rename-reads_complete; then
  # Delete first two read files because they are technical
  rm -f ${sra_id}_{1,2}.fastq
  # Rename the read file with biological reads
  mv ./${sra_id}_3.fastq ./${sra_id}_S1_L001_R1_001.fastq 
  # Rename the read file with index
  mv ./${sra_id}_4.fastq ./${sra_id}_S1_L001_R2_001.fastq
  
  touch ./${sra_id}_status/3_rename-reads_complete
fi

# Perform sequence alignment
if ! test -f ./${sra_id}_status/4_alignment_complete; then
  # Delete incomplete files if they exist
  if test -f ./${sra_id}-ensembl; then
    rm -r ./${sra_id}-ensembl
  fi
  # Align sequences with cell ranger
  cellranger count --id=${sra_id}-ensembl \
                 --transcriptome=${transcriptome_ref_path} \
                 --fastqs=. \
                 --sample=${sra_id} \
                 --localcores=10 \
                 --localmem=100
                 
  touch ./${sra_id}_status/4_alignment_complete
fi


# Save alignment output to parent directory
if ! test -f ./${sra_id}_status/5_save-alignment_complete; then
  # Copy the filtered matrix, zip and delete
  mv ./${sra_id}-ensembl/outs/filtered_feature_bc_matrix ../${sra_id}-ensembl_filtered_feature_bc_matrix
  tar -czvf ../${sra_id}-ensembl_filtered_feature_bc_matrix.tar.gz ../${sra_id}-ensembl_filtered_feature_bc_matrix
  rm -r ../${sra_id}-ensembl_filtered_feature_bc_matrix
  
  touch ./${sra_id}_status/5_save-alignment_complete
fi


# Delete downloaded and intermediate data
if ! test -f ./${sra_id}_status/6_cleanup_complete; then
  # Delete reads, cellranger output, and SRA files
  rm -f ./${sra_id}_S1_L001_R{1,2}_001.fastq  
  rm -r ./${sra_id}-ensembl
  rm -r ./${sra_id}
  touch ./${sra_id}_status/6_cleanup_complete
fi

rm -r ./${sra_id}_status





