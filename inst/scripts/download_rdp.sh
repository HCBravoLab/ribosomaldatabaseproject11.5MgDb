## Download bacterial and archaea unaligned seuqences
wget --no-check-certificate https://rdp.cme.msu.edu/download/current_Bacteria_unaligned.fa.gz
wget --no-check-certificate https://rdp.cme.msu.edu/download/current_Archaea_unaligned.fa.gz

cat current_Bacteria_unaligned.fa.gz current_Archaea_unaligned.fa.gz > current_Prokaryota_unaligned.fa.gz
#split -b 45m current_Bacteria_unaligned.fa.gz current_Bacteria_unaligned_subset
