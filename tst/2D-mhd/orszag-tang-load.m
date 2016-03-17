%this works with change to line 14-15 in ath_parse_filename
%[path, file, ext, versn] = fileparts(filename);
%[path, file, ext] = fileparts(filename);


[mygrid,status]=ath_init_grid('OrszagTang.0071.bin');
[time,dt,rho,status]=ath_readbin(mygrid,'OrszagTang.0071.bin','d');
surf(rho);
surf(rho,'LineStyle','none');