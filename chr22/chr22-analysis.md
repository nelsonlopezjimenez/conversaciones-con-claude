# Prompt
```sh
(genomicschr22) usergolden26ubuntu2@CISGOLDEN26:~/chr22$ ls
ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
gencode.v19.annotation.gtf
human_g1k_v37.fasta
integrated_call_samples_v3.20130502.ALL.panel; python script combination ok, i need to install some dependencies though; step-by-step preferred; jypyter notebook maybe later;
```

Perfect! You have everything needed. Let's start step-by-step.

## Step 0: Install Dependencies

```bash
# Activate your environment
conda activate genomicschr22

# Install required packages
conda install -c bioconda pysam bcftools samtools
conda install -c conda-forge biopython pandas

# Verify installations
python -c "import pysam; import Bio; import pandas; print('All packages installed successfully')"
```

---

## Step 1: Parse GTF and Extract Protein-Coding Exons

**Purpose:** Extract all exon coordinates for protein-coding genes on chromosome 22

**Create file:** `step1_parse_gtf.py`

```python
#!/usr/bin/env python3
"""
Step 1: Parse GTF and extract protein-coding exons from chromosome 22
"""
import pandas as pd
import json
from collections import defaultdict

def parse_gtf_exons(gtf_file, chromosome='22', output_json='chr22_exons.json'):
    """
    Extract protein-coding exons from GTF file
    
    Returns dictionary: {gene_id: {gene_name, transcript_id, exons: [(start, end, exon_number)]}}
    """
    print(f"Parsing {gtf_file}...")
    
    genes_data = defaultdict(lambda: {
        'gene_name': '',
        'transcripts': defaultdict(list)
    })
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom = fields[0]
            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            attributes = fields[8]
            
            # Only chromosome 22
            if chrom != chromosome and chrom != f'chr{chromosome}':
                continue
            
            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr:
                    try:
                        key, value = attr.split(' ', 1)
                        attr_dict[key] = value.strip('"')
                    except:
                        pass
            
            # Only protein-coding genes
            if attr_dict.get('gene_type') != 'protein_coding':
                continue
            
            gene_id = attr_dict.get('gene_id', '')
            gene_name = attr_dict.get('gene_name', '')
            
            if feature == 'exon':
                transcript_id = attr_dict.get('transcript_id', '')
                exon_number = attr_dict.get('exon_number', '0')
                
                genes_data[gene_id]['gene_name'] = gene_name
                genes_data[gene_id]['transcripts'][transcript_id].append({
                    'start': start,
                    'end': end,
                    'exon_number': int(exon_number)
                })
    
    # Select canonical transcript (longest) for each gene
    canonical_genes = {}
    
    for gene_id, data in genes_data.items():
        gene_name = data['gene_name']
        
        # Find longest transcript
        longest_transcript = None
        max_length = 0
        
        for transcript_id, exons in data['transcripts'].items():
            total_length = sum(e['end'] - e['start'] + 1 for e in exons)
            if total_length > max_length:
                max_length = total_length
                longest_transcript = transcript_id
        
        if longest_transcript:
            exons = sorted(data['transcripts'][longest_transcript], 
                          key=lambda x: x['exon_number'])
            
            canonical_genes[gene_id] = {
                'gene_name': gene_name,
                'transcript_id': longest_transcript,
                'exon_count': len(exons),
                'exons': exons
            }
    
    print(f"\nFound {len(canonical_genes)} protein-coding genes")
    print(f"Total exons: {sum(g['exon_count'] for g in canonical_genes.values())}")
    
    # Save to JSON
    with open(output_json, 'w') as f:
        json.dump(canonical_genes, f, indent=2)
    
    print(f"Saved to {output_json}")
    
    # Print summary statistics
    exon_counts = [g['exon_count'] for g in canonical_genes.values()]
    print(f"\nExon count statistics:")
    print(f"  Min: {min(exon_counts)}")
    print(f"  Max: {max(exon_counts)}")
    print(f"  Median: {sorted(exon_counts)[len(exon_counts)//2]}")
    
    return canonical_genes

if __name__ == '__main__':
    genes = parse_gtf_exons('gencode.v19.annotation.gtf')
    
    # Show first 5 genes as example
    print("\nExample genes:")
    for i, (gene_id, data) in enumerate(list(genes.items())[:5]):
        print(f"\n{data['gene_name']} ({gene_id}):")
        print(f"  Transcript: {data['transcript_id']}")
        print(f"  Exons: {data['exon_count']}")
        print(f"  First exon: {data['exons'][0]['start']}-{data['exons'][0]['end']}")
```

**Run it:**
```bash
python step1_parse_gtf.py
```

**Expected output:**
- `chr22_exons.json` file created
- Statistics about genes and exons found
- Should show ~500-600 protein-coding genes
- ~4,500-5,000 total exons

---

**Does this run successfully?** Once confirmed, I'll provide Step 2 (extracting trio genotypes from VCF).