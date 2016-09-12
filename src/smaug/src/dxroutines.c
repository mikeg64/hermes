#include <stdio.h>

void gendxgen(char *dir,char *jobname,int nsteps,int n1,int n2)
{
//generate dx general file 
//describing magnetic fields
//particle locations
  //file=out/jobform.out
  //grid = 1
  //format = ascii
  //interleaving = record
  //majority = row
  //field = nsteps, nx, ny
  //structure = scalar, scalar, scalar
  //type = int, int, int
  //dependency = positions, positions, positions
  //positions = regular, 0, 1

  //end
  char dxformgenfile[300];
  char dxgenfile[300];
  char basedxgenfile[300];
//printf("here in gendxgen1 \n"); 
  sprintf(dxformgenfile,"dx/%s_form.general",jobname);
//printf("here in gendxgen1a %s %s\n",dxformgenfile,dir); 
  FILE *fdform=fopen(dxformgenfile,"w");
    fprintf(fdform, "file=%s/form%s.out\n",dir,jobname);
    fprintf(fdform,"grid=1\n");
    fprintf(fdform,"format = ascii \n interleaving = record \n majority = row \n");
    fprintf(fdform, "field = nsteps, nx, ny \n structure = scalar, scalar, scalar \n type = int, int, int  \n dependency = positions, positions,positions  \n positions = regular, 0, 1 \n end \n ");
  fclose(fdform);    
//generate dx general file for this data set
  //file=out/job.out
  //grid 51 x 51
  //format = ascii
  //interleaving = field
  //majority = row
  //header = lines 1

  //series =  24 , 1, 1, separator=lines 1
  //field = field0, field1
  //structure = 2-vector, scalar
  //type = float, float
  //dependency = positions, positions
  //positions = regular,regular, 0, 1,0,1

  //end
printf("here in gendxgen2 \n");
  sprintf(dxgenfile,"dx/%s.general",jobname);

  fdform=fopen(dxgenfile,"w");
    fprintf(fdform, "file=%s/%s.out\n",dir,jobname);
    fprintf(fdform,"grid %d X %d\n",n1,n2);
    fprintf(fdform,"format = ascii \n interleaving = field \n majority = row \n header = lines 1 \n");
    fprintf(fdform, "series =  %d  , 1, 1, separator=lines 1\n",nsteps-1);
    fprintf(fdform, "field = field0, field1 \n structure = 2-vector, scalar \n type = float, float  \n dependency = positions, positions  \n positions = regular,regular, 0, 1,0,1 \n end \n ");
  fclose(fdform);


 sprintf(basedxgenfile,"dx/base%s.general",jobname);

  FILE *fdbform=fopen(basedxgenfile,"w");
   // mfprintf(fdbform, 'file=%s\n', directory+'/'+jobname+'.out');
    fprintf(fdbform,"\ngrid %d X %d\n",n1,n2);
    fprintf(fdbform,"format = ascii \n interleaving = field \n majority = row \n header = lines 1 \n");
    fprintf(fdbform, "series =  1  , 1, 1, separator=lines 1\n");
    fprintf(fdbform, "field = field0, field1 \n structure = 2-vector, scalar \n type = float, float  \n dependency = positions, positions  \n positions = regular,regular, 0, 1,0,1 \n end \n ");
  fclose(fdbform);

printf("here in gendxgen3 \n");

//endfunction
}
