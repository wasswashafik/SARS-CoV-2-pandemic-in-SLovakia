
[A,infectedString] = system('curl http://covid.gergelytakacs.com/csv.php?dataset=infected -s');
fid = fopen('infected.csv','wt');
fprintf(fid, infectedString);
fclose(fid);
clear A infectedString fid

[A,recoveredString] = system('curl http://covid.gergelytakacs.com/csv.php?dataset=recovered -s');
fid = fopen('recovered.csv','wt');
fprintf(fid, recoveredString);
fclose(fid);
clear A recoveredString fid

[A,deceasedString] = system('curl http://covid.gergelytakacs.com/csv.php?dataset=deceased -s');
fid = fopen('deceased.csv','wt');
fprintf(fid, deceasedString);
fclose(fid);
clear A deceasedString fid
