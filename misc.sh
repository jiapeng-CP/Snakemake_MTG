
# merge tables from individual samples for MetaPhlan
~/bin/biobakery4/bin/merge_metaphlan_tables.py MetaPhlan/*txt > merged.metaphlan.txt



# merge tables from individual samples for HumanN
/home/artemisl/.conda/envs/biobakery/bin/humann_join_tables --input HumanN --search-subdirectories --file_name genefamilies.tsv --output merge.humann.genef.tsv &
/home/artemisl/.conda/envs/biobakery/bin/humann_join_tables --input HumanN --search-subdirectories --file_name pathabundance.tsv --output merge.humann.patha.tsv &
/home/artemisl/.conda/envs/biobakery/bin/humann_join_tables --input HumanN --search-subdirectories --file_name pathcoverage.tsv --output merge.humann.pathc.tsv &



# produce individual funtional tables
humann_regroup_table -i merge.humann.genef.tsv -g uniref90_ko -o humann.ko.tsv &
humann_regroup_table -i merge.humann.genef.tsv -g uniref90_rxn -o humann.rxn.tsv &
humann_regroup_table -i merge.humann.genef.tsv -g uniref90_go -o humann.go.tsv &
humann_regroup_table -i merge.humann.genef.tsv -g uniref90_level4ec -o humann.level4ec.tsv &
humann_regroup_table -i merge.humann.genef.tsv -g uniref90_pfam -o humann.pfam.tsv &
humann_regroup_table -i merge.humann.genef.tsv -g uniref90_eggnog -o humann.eggnog.tsv &
