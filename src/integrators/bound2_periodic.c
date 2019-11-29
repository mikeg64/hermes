//boundary.c

 
k=0;    
for(dim=0; dim<2; dim++)
{
        for (ib=0; ib<=1; ib++)   //go through upper boundary and lower boundary
    {
            
            
            
            
            
            
            
      switch(dim)
      {
          case 0:
             if(ib==0)
                 //for(i1=0;i1<n1z;i1++)
                 //{
                    for(j=0;j<n2z;j++)
                    {
                            //for(i3=0;i1<n3z;i3++)
                            //{

        Uinit[k][j][n1z-1].d  = Uinit[k][j][3].d;
        Uinit[k][j][n1z-1].M1 = Uinit[k][j][3].M1;
        Uinit[k][j][n1z-1].M2 = Uinit[k][j][3].M2;
        Uinit[k][j][n1z-1].M3 = Uinit[k][j][3].M3;
#ifndef BAROTROPIC
        Uinit[k][j][n1z-1].E  = Uinit[k][j][3].E;
#endif /* BAROTROPIC */
#ifdef MHD
        Uinit[k][j][n1z-1].B1c = Uinit[k][j][3].B1c;
        Uinit[k][j][n1z-1].B2c = Uinit[k][j][3].B2c;
        Uinit[k][j][n1z-1].B3c = Uinit[k][j][3].B3c;
        //Bxc[i] = Uinit[k][j][i].B1c;
        //Bxb[i] = pG->B1cb[k][j][i];
        //B1_x1Face[k][j][i] = pG->B1i[k][j][i];
#endif /* MHD */


#ifdef SAC_INTEGRATOR
        Uinit[k][j][n1z-1].db  = Uinit[k][j][3].db;
#ifdef MHD
        Uinit[k][j][n1z-1].B1cb = Uinit[k][j][3].B1cb;
        Uinit[k][j][n1z-1].B2cb = Uinit[k][j][3].B2cb;
        Uinit[k][j][n1z-1].B3cb = Uinit[k][j][3].B3cb;
#endif

#endif
#ifdef SMAUG_INTEGRATOR
        Uinit[k][j][n1z-1].db  = Uinit[k][j][3].db;

#ifdef MHD
        Uinit[k][j][n1z-1].B1cb = Uinit[k][j][3].B1cb;
        Uinit[k][j][n1z-1].B2cb = Uinit[k][j][3].B2cb;
        Uinit[k][j][n1z-1].B3cb = Uinit[k][j][3].B3cb;
#endif 
#endif                        
                        
                    
                    
                    
                    
                    
        Uinit[k][j][n1z-2].d  = Uinit[k][j][2].d;
        Uinit[k][j][n1z-2].M1 = Uinit[k][j][2].M1;
        Uinit[k][j][n1z-2].M2 = Uinit[k][j][2].M2;
        Uinit[k][j][n1z-2].M3 = Uinit[k][j][2].M3;
#ifndef BAROTROPIC
        Uinit[k][j][n1z-2].E  = Uinit[k][j][2].E;
#endif /* BAROTROPIC */
#ifdef MHD
        Uinit[k][j][n1z-2].B1c = Uinit[k][j][2].B1c;
        Uinit[k][j][n1z-2].B2c = Uinit[k][j][2].B2c;
        Uinit[k][j][n1z-2].B3c = Uinit[k][j][2].B3c;
        //Bxc[i] = Uinit[k][j][i].B1c;
        //Bxb[i] = pG->B1cb[k][j][i];
        //B1_x1Face[k][j][i] = pG->B1i[k][j][i];
#endif /* MHD */


#ifdef SAC_INTEGRATOR
        Uinit[k][j][n1z-2].db  = Uinit[k][j][2].db;
#ifdef MHD
        Uinit[k][j][n1z-2].B1cb = Uinit[k][j][2].B1cb;
        Uinit[k][j][n1z-2].B2cb = Uinit[k][j][2].B2cb;
        Uinit[k][j][n1z-2].B3cb = Uinit[k][j][2].B3cb;
#endif

#endif
#ifdef SMAUG_INTEGRATOR
        Uinit[k][j][n1z-2].db  = Uinit[k][j][2].db;

#ifdef MHD
        Uinit[k][j][n1z-2].B1cb = Uinit[k][j][2].B1cb;
        Uinit[k][j][n1z-2].B2cb = Uinit[k][j][2].B2cb;
        Uinit[k][j][n1z-2].B3cb = Uinit[k][j][2].B3cb;
#endif
#endif                  
                    
 
                                
                            //}
                    }
                 //}
      
    
                    
                                   
                    
             else
                 //for(i1=0;i1<n1z;i1++)
                 //{
                    for(j=0;j<n2z;j++)
                    {
                            //for(i3=0;i1<n3z;i3++)
                            //{


        Uinit[k][j][0].d  = Uinit[k][j][n1z-4].d;
        Uinit[k][j][0].M1 = Uinit[k][j][n1z-4].M1;
        Uinit[k][j][0].M2 = Uinit[k][j][n1z-4].M2;
        Uinit[k][j][0].M3 = Uinit[k][j][n1z-4].M3;
#ifndef BAROTROPIC
        Uinit[k][j][0].E  = Uinit[k][j][n1z-4].E;
#endif /* BAROTROPIC */
#ifdef MHD
        Uinit[k][j][0].B1c = Uinit[k][j][n1z-4].B1c;
        Uinit[k][j][0].B2c = Uinit[k][j][n1z-4].B2c;
        Uinit[k][j][0].B3c = Uinit[k][j][n1z-4].B3c;
        //Bxc[i] = Uinit[k][j][i].B1c;
        //Bxb[i] = pG->B1cb[k][j][i];
        //B1_x1Face[k][j][i] = pG->B1i[k][j][i];
#endif /* MHD */


#ifdef SAC_INTEGRATOR
        Uinit[k][j][0].db  = Uinit[k][j][n1z-4].db;
#ifdef MHD
        Uinit[k][j][0].B1cb = Uinit[k][j][n1z-4].B1cb;
        Uinit[k][j][0].B2cb = Uinit[k][j][n1z-4].B2cb;
        Uinit[k][j][0].B3cb = Uinit[k][j][n1z-4].B3cb;
#endif

#endif
#ifdef SMAUG_INTEGRATOR
        Uinit[k][j][0].db  = Uinit[k][j][n1z-4].db;

#ifdef MHD
        Uinit[k][j][0].B1cb = Uinit[k][j][n1z-4].B1cb;
        Uinit[k][j][0].B2cb = Uinit[k][j][n1z-4].B2cb;
        Uinit[k][j][0].B3cb = Uinit[k][j][n1z-4].B3cb;
#endif 
#endif



        Uinit[k][j][1].d  = Uinit[k][j][n1z-3].d;
        Uinit[k][j][1].M1 = Uinit[k][j][n1z-3].M1;
        Uinit[k][j][1].M2 = Uinit[k][j][n1z-3].M2;
        Uinit[k][j][1].M3 = Uinit[k][j][n1z-3].M3;
#ifndef BAROTROPIC
        Uinit[k][j][1].E  = Uinit[k][j][n1z-3].E;
#endif /* BAROTROPIC */
#ifdef MHD
        Uinit[k][j][1].B1c = Uinit[k][j][n1z-3].B1c;
        Uinit[k][j][1].B2c = Uinit[k][j][n1z-3].B2c;
        Uinit[k][j][1].B3c = Uinit[k][j][n1z-3].B3c;
        //Bxc[i] = Uinit[k][j][i].B1c;
        //Bxb[i] = pG->B1cb[k][j][i];
        //B1_x1Face[k][j][i] = pG->B1i[k][j][i];
#endif /* MHD */


#ifdef SAC_INTEGRATOR
        Uinit[k][j][1].db  = Uinit[k][j][n1z-3].db;
#ifdef MHD
        Uinit[k][j][1].B1cb = Uinit[k][j][n1z-3].B1cb;
        Uinit[k][j][1].B2cb = Uinit[k][j][n1z-3].B2cb;
        Uinit[k][j][1].B3cb = Uinit[k][j][n1z-3].B3cb;
#endif

#endif
#ifdef SMAUG_INTEGRATOR
        Uinit[k][j][1].db  = Uinit[k][j][n1z-3].db;

#ifdef MHD
        Uinit[k][j][1].B1cb = Uinit[k][j][n1z-3].B1cb;
        Uinit[k][j][1].B2cb = Uinit[k][j][n1z-3].B2cb;
        Uinit[k][j][1].B3cb = Uinit[k][j][n1z-3].B3cb;
#endif 
#endif



                                
                            //}
                    }
                 //}
              
                   
              
          case 1:
             if(ib==0)
                 for(i=0;i<n1z;i++)
                 {
                    //for(i2=0;i2<n2z;i2++)
                    //{
                            //for(i3=0;i1<n3z;i3++)
                            //{




        Uinit[k][n2z-1][i].d  = Uinit[k][3][i].d;
        Uinit[k][n2z-1][i].M1 = Uinit[k][3][i].M1;
        Uinit[k][n2z-1][i].M2 = Uinit[k][3][i].M2;
        Uinit[k][n2z-1][i].M3 = Uinit[k][3][i].M3;
#ifndef BAROTROPIC
        Uinit[k][n2z-1][i].E  = Uinit[k][3][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        Uinit[k][n2z-1][i].B1c = Uinit[k][3][i].B1c;
        Uinit[k][n2z-1][i].B2c = Uinit[k][3][i].B2c;
        Uinit[k][n2z-1][i].B3c = Uinit[k][3][i].B3c;
        //Bxc[i] = Uinit[k][j][i].B1c;
        //Bxb[i] = pG->B1cb[k][j][i];
        //B1_x1Face[k][j][i] = pG->B1i[k][j][i];
#endif /* MHD */


#ifdef SAC_INTEGRATOR
        Uinit[k][n2z-1][i].db  = Uinit[k][3][i].db;
#ifdef MHD
        Uinit[k][n2z-1][i].B1cb = Uinit[k][3][i].B1cb;
        Uinit[k][n2z-1][i].B2cb = Uinit[k][3][i].B2cb;
        Uinit[k][n2z-1][i].B3cb = Uinit[k][3][i].B3cb;
#endif

#endif
#ifdef SMAUG_INTEGRATOR
        Uinit[k][n2z-1][i].db  = Uinit[k][3][i].db;

#ifdef MHD
        Uinit[k][n2z-1][i].B1cb = Uinit[k][3][i].B1cb;
        Uinit[k][n2z-1][i].B2cb = Uinit[k][3][i].B2cb;
        Uinit[k][n2z-1][i].B3cb = Uinit[k][3][i].B3cb;
#endif 
#endif                    
             
                    
                    
        Uinit[k][n2z-2][i].d  =  Uinit[k][2][i].d;
        Uinit[k][n2z-2][i].M1 =  Uinit[k][2][i].M1;
        Uinit[k][n2z-2][i].M2 =  Uinit[k][2][i].M2;
        Uinit[k][n2z-2][i].M3 =  Uinit[k][2][i].M3;
#ifndef BAROTROPIC
        Uinit[k][n2z-2][i].E  =  Uinit[k][2][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        Uinit[k][n2z-2][i].B1c =  Uinit[k][2][i].B1c;
        Uinit[k][n2z-2][i].B2c =  Uinit[k][2][i].B2c;
        Uinit[k][n2z-2][i].B3c =  Uinit[k][2][i].B3c;
        //Bxc[i] = Uinit[k][j][i].B1c;
        //Bxb[i] = pG->B1cb[k][j][i];
        //B1_x1Face[k][j][i] = pG->B1i[k][j][i];
#endif /* MHD */


#ifdef SAC_INTEGRATOR
        Uinit[k][n2z-2][i].db  =  Uinit[k][2][i].db;
#ifdef MHD
        Uinit[k][n2z-2][i].B1cb =  Uinit[k][2][i].B1cb;
        Uinit[k][n2z-2][i].B2cb =  Uinit[k][2][i].B2cb;
        Uinit[k][n2z-2][i].B3cb =  Uinit[k][2][i].B3cb;
#endif

#endif
#ifdef SMAUG_INTEGRATOR
        Uinit[k][n2z-2][i].db  =  Uinit[k][2][i].db;

#ifdef MHD
        Uinit[k][n2z-2][i].B1cb =  Uinit[k][2][i].B1cb;
        Uinit[k][n2z-2][i].B2cb =  Uinit[k][2][i].B2cb;
        Uinit[k][n2z-2][i].B3cb =  Uinit[k][2][i].B3cb;
#endif                             
#endif
                          
                                
                            //}
                    //}
                 }
                 
             else
                 for(i=0;i<n1z;i++)
                 {
                    //for(i2=0;i2<n2z;i2++)
                    //{
                            //for(i3=0;i1<n3z;i3++)
                            //{




         Uinit[k][0][i].d  =  Uinit[k][n2z-4][i].d;
         Uinit[k][0][i].M1 =  Uinit[k][n2z-4][i].M1;
         Uinit[k][0][i].M2 =  Uinit[k][n2z-4][i].M2;
         Uinit[k][0][i].M3 =  Uinit[k][n2z-4][i].M3;
#ifndef BAROTROPIC
         Uinit[k][0][i].E  =  Uinit[k][n2z-4][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
         Uinit[k][0][i].B1c =  Uinit[k][n2z-4][i].B1c;
         Uinit[k][0][i].B2c =  Uinit[k][n2z-4][i].B2c;
         Uinit[k][0][i].B3c =  Uinit[k][n2z-4][i].B3c;
        //Bxc[i] = Uinit[k][j][i].B1c;
        //Bxb[i] = pG->B1cb[k][j][i];
        //B1_x1Face[k][j][i] = pG->B1i[k][j][i];
#endif /* MHD */


#ifdef SAC_INTEGRATOR
         Uinit[k][0][i].db  =  Uinit[k][n2z-4][i].db;
#ifdef MHD
         Uinit[k][0][i].B1cb =  Uinit[k][n2z-4][i].B1cb;
         Uinit[k][0][i].B2cb =  Uinit[k][n2z-4][i].B2cb;
         Uinit[k][0][i].B3cb =  Uinit[k][n2z-4][i].B3cb;
#endif
#endif
                     
#ifdef SMAUG_INTEGRATOR
         Uinit[k][0][i].db  =  Uinit[k][n2z-4][i].db;

#ifdef MHD
         Uinit[k][0][i].B1cb =  Uinit[k][n2z-4][i].B1cb;
         Uinit[k][0][i].B2cb =  Uinit[k][n2z-4][i].B2cb;
         Uinit[k][0][i].B3cb =  Uinit[k][n2z-4][i].B3cb;
#endif 
#endif



         Uinit[k][1][i].d  =  Uinit[k][n2z-3][i].d;
         Uinit[k][1][i].M1 =  Uinit[k][n2z-3][i].M1;
         Uinit[k][1][i].M2 =  Uinit[k][n2z-3][i].M2;
         Uinit[k][1][i].M3 =  Uinit[k][n2z-3][i].M3;
#ifndef BAROTROPIC
         Uinit[k][1][i].E  =  Uinit[k][n2z-3][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
         Uinit[k][1][i].B1c =  Uinit[k][n2z-3][i].B1c;
         Uinit[k][1][i].B2c =  Uinit[k][n2z-3][i].B2c;
         Uinit[k][1][i].B3c =  Uinit[k][n2z-3][i].B3c;
        //Bxc[i] = Uinit[k][j][i].B1c;
        //Bxb[i] = pG->B1cb[k][j][i];
        //B1_x1Face[k][j][i] = pG->B1i[k][j][i];
#endif /* MHD */


#ifdef SAC_INTEGRATOR
         Uinit[k][1][i].db  =  Uinit[k][n2z-3][i].db;
#ifdef MHD
         Uinit[k][1][i].B1cb =  Uinit[k][n2z-3][i].B1cb;
         Uinit[k][1][i].B2cb =  Uinit[k][n2z-3][i].B2cb;
         Uinit[k][1][i].B3cb =  Uinit[k][n2z-3][i].B3cb;
#endif
#endif
                     
#ifdef SMAUG_INTEGRATOR
         Uinit[k][1][i].db  =  Uinit[k][n2z-3][i].db;

#ifdef MHD
         Uinit[k][1][i].B1cb =  Uinit[k][n2z-3][i].B1cb;
         Uinit[k][1][i].B2cb =  Uinit[k][n2z-3][i].B2cb;
         Uinit[k][1][i].B3cb =  Uinit[k][n2z-3][i].B3cb;
#endif 
#endif
                              
                                
                            //}
                    //}
                 }
         
                                        
              
              
      }   //end of switch(dim)
        
        
        
    }    //end of ib loop
} //looping over dimensions    
    
    
    
    
    
    
    
    
    
    

