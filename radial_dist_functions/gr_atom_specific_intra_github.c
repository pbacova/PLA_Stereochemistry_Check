//Intramolecular radial distribution functions for specific atoms in PLA polymer. Check 100LA.itp file in stereochemistry_check directory for the order of the atoms in the gro file.
//If the input file contain also velocities, uncomment line 157 and comment the line 156.
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include<string.h>

int nchain=13;
int nframes=2401;
int lchain=321;
int am=9;
double grid, boxmaxx, boxmaxy,boxmaxz;
int skip=0;
char **p;
void print_help()
{
	printf ("usage:re [options]\n \
	options:\n \
	-f number of frames\n \
	-p number of chains\n \
        -k number of atoms in every chain\n \
        -o number of frames to skip\n \
        -g grid(default 0.1)\n \
        -x boxmax x\n \
        -y boxmax y\n \
        -z boxmax z\n \
	-h help\n");
	exit(1);
}

int parse_options (int argc, char *argv[])
{
	if (argc==1) print_help();
	while (--argc && (*++argv)[0] == '-') 
	{
		char c;
		c = *++argv[0];
		switch (c) 
		{
			case 'f':
				nframes=atoi(*++argv);
				break;
			case 'p':
				nchain=atoi(*++argv);
				break;
			case 'k':
				lchain=atoi(*++argv);
				break;
                        case 'o':
				skip=atoi(*++argv);
				break;
			case 'g':
				grid=strtod(*++argv,p);
				break;
			case 'x':
				boxmaxx=strtod(*++argv,p);
				break;
			case 'y':
				boxmaxy=strtod(*++argv,p);
				break;
			case 'z':
				boxmaxz=strtod(*++argv,p);
				break;
			case 'h' :
				print_help();
				break;
			default:
				printf ("\ncorre: incorrect input, help:-h \
					\n\n", argv[0]);
				exit (1);
		
		}
		--argc;
	}
	return 0;
}
int main(int argc, char *argv[]){
	// initialize
	double x,y,z,vx,vy,vz;
	FILE *init,*out1;
	int  i1,i2,ir,o,i,j,k,a,b,c,*nmol,*natom,start,end,start2,end2;
        int ngridwr;
        char outname[150];
        char buf[100]; 
	char **string_array;
        double *xm,*ym,*zm,rho;
        double ifac,rx,ry,rz,rr,*grO1,*grC2,*grC3,*grO4,*grH5,*grC6,*grH7,*grH10,*grO10,*grH11,maxd,box2,boxx,boxy,boxz;
        char tmp1[5];
        parse_options (argc, argv);
	
         if((init=fopen("input.gro","r"))==NULL){
           printf("file can not be opened\n");
           return 1;
           } 
          int iframes=nframes-skip; 
          int N_all=nchain*lchain;
          int nmon=lchain/am-2;
          box2 = boxmaxx/2;
          maxd = box2*sqrt(3);
	  ngridwr = (int)(maxd/grid);
	  long int numO1[ngridwr];
	  long int numC2[ngridwr];
	  long int numC3[ngridwr];
	  long int numO4[ngridwr];
	  long int numH5[ngridwr];
	  long int numC6[ngridwr];
	  long int numH7[ngridwr];
	  long int numO10[ngridwr];
	  long int numH10[ngridwr];
	  long int numH11[ngridwr];
  	  grO1=(double*)malloc((ngridwr)*sizeof(double));
  	  grC2=(double*)malloc((ngridwr)*sizeof(double));
  	  grC3=(double*)malloc((ngridwr)*sizeof(double));
  	  grO4=(double*)malloc((ngridwr)*sizeof(double));
  	  grH5=(double*)malloc((ngridwr)*sizeof(double));
  	  grC6=(double*)malloc((ngridwr)*sizeof(double));
  	  grH7=(double*)malloc((ngridwr)*sizeof(double));
  	  grO10=(double*)malloc((ngridwr)*sizeof(double));
  	  grH10=(double*)malloc((ngridwr)*sizeof(double));
  	  grH11=(double*)malloc((ngridwr)*sizeof(double));
          for (i=0;i<ngridwr;i++) grO1[i]=grC2[i]=grC3[i]=grO4[i]=grH5[i]=grC6[i]=grH7[i]=grO10[i]=grH10[i]=grH11[i]=0.0;
          rho=(4.0/3.0)*M_PI*N_all*((grid*grid*grid)/(boxmaxx*boxmaxy*boxmaxz));
         
           sprintf(outname,"gr_atom_specific_intra");
          out1=fopen(outname,"w");
          xm=(double*)malloc((N_all)*sizeof(double));
          ym=(double*)malloc((N_all)*sizeof(double));
          zm=(double*)malloc((N_all)*sizeof(double));
          nmol=(int*)malloc(N_all*sizeof(int));
          natom=(int*)malloc(N_all*sizeof(int));
          string_array=(char**)malloc(N_all*sizeof(char*));
         for(i=0;i<N_all;i++){
          string_array[i]=(char*)malloc(5*sizeof(char));
          nmol[i]=natom[i]=0;
          }
         for(j=0;j<N_all;j++)xm[j]=ym[j]=zm[j]=0.0; 
for(j=0;j<skip;j++){
        fgets(buf,sizeof(buf),init);
        fgets(buf,sizeof(buf),init);
        for(i=0;i<N_all;i++){
         fgets(buf,sizeof(buf),init);
         }
        fscanf(init,"%lf%lf%lf",&boxx,&boxy,&boxz);
        fgets(buf,sizeof(buf),init);
}
 for(j=0;j<(iframes);j++){
	for (i=0; i<ngridwr; i++) numO1[i]=numC2[i]=numC3[i]=numO4[i]=numH5[i]=numC6[i]=numH7[i]=numO10[i]=numH10[i]=numH11[i]=0;
        fgets(buf,sizeof(buf),init);
        fgets(buf,sizeof(buf),init);
    for(a=0;a<nchain;a++){ 
        //chain 
           start=a*(lchain) ;
           end=start+lchain;
           for(i=start;i<end;i++){
             fscanf(init,"%5d%5c%5c%5d%lf%lf%lf",&nmol[i],&tmp1[0],string_array[i],&natom[i],&xm[i],&ym[i],&zm[i]);
            // fscanf(init,"%5d%5c%5c%5d%lf%lf%lf%lf%lf%lf",&nmol[i],&tmp1[0],string_array[i],&natom[i],&xm[i],&ym[i],&zm[i],&vx,&vy,&vz);
             }
       }//end of reading chains
       fscanf(init,"%lf%lf%lf",&boxx,&boxy,&boxz);
       fgets(buf,sizeof(buf),init);
       for (i=0; i<N_all; i++){
		xm[i] = xm[i] - boxx*((int)(xm[i]/boxx));
		ym[i] = ym[i] - boxy*((int)(ym[i]/boxy));
		zm[i] = zm[i] - boxz*((int)(zm[i]/boxz));
		if (xm[i]<0) xm[i] = xm[i] + boxx;
		if (ym[i]<0) ym[i] = ym[i] + boxy;
		if (zm[i]<0) zm[i] = zm[i] + boxz;
	}
    for(a=0;a<(nchain-1);a++){ 
    //check the first monomer in the chain
          start=a*lchain;
          end=start+10;
          i2=end-3;
          //intramonomer interactions: loop over hydrogens
          for(i=i2;i<end;i++){
          //O1
            i1=start;
            //printf("%s %s\n",string_array[i1],string_array[i]);
            //if(a==nchain-1)printf("%lf %lf %lf\n",xm[i],ym[i],zm[i]);
   	    rx = ( xm[i1] - xm[i]) - boxx* ((int)(2*(xm[i1]-xm[i])/boxx));
 	    ry = ( ym[i1] - ym[i]) - boxy* ((int)(2*(ym[i1]-ym[i])/boxy));
 	    rz = ( zm[i1] - zm[i]) - boxz* ((int)(2*(zm[i1]-zm[i])/boxz));
 	    rr = sqrt(rx*rx+ry*ry+rz*rz);
	    ir = (int)(rr/grid);
//            printf("%s\n",string_array[i1]);
	    numO1[ir]=numO1[ir]+1;
          //O4
            i1=start+3;
   	    rx = ( xm[i1] - xm[i]) - boxx* ((int)(2*(xm[i1]-xm[i])/boxx));
 	    ry = ( ym[i1] - ym[i]) - boxy* ((int)(2*(ym[i1]-ym[i])/boxy));
 	    rz = ( zm[i1] - zm[i]) - boxz* ((int)(2*(zm[i1]-zm[i])/boxz));
 	    rr = sqrt(rx*rx+ry*ry+rz*rz);
	    ir = (int)(rr/grid);
            //printf("%s\n",string_array[i1]);
	    numO4[ir]=numO4[ir]+1;
          //H10
            i1=start+5;
   	    rx = ( xm[i1] - xm[i]) - boxx* ((int)(2*(xm[i1]-xm[i])/boxx));
 	    ry = ( ym[i1] - ym[i]) - boxy* ((int)(2*(ym[i1]-ym[i])/boxy));
 	    rz = ( zm[i1] - zm[i]) - boxz* ((int)(2*(zm[i1]-zm[i])/boxz));
 	    rr = sqrt(rx*rx+ry*ry+rz*rz);
	    ir = (int)(rr/grid);
            //printf("%s\n",string_array[i1]);
	    numH10[ir]=numH10[ir]+1;
           }
          //H10 H5; H10 C6   
           i1=start+5;
           i2=start+4;
   	   rx = ( xm[i1] - xm[i2]) - boxx* ((int)(2*(xm[i1]-xm[i2])/boxx));
 	   ry = ( ym[i1] - ym[i2]) - boxy* ((int)(2*(ym[i1]-ym[i2])/boxy));
 	   rz = ( zm[i1] - zm[i2]) - boxz* ((int)(2*(zm[i1]-zm[i2])/boxz));
 	   rr = sqrt(rx*rx+ry*ry+rz*rz);
	   ir = (int)(rr/grid);
          // printf("%s\n",string_array[i1]);
	   numH10[ir]=numH10[ir]+1;
           i1=start+5;
           i2=start+6;
           // printf("%s %s\n",string_array[i1],string_array[i2]);
   	   rx = ( xm[i1] - xm[i2]) - boxx* ((int)(2*(xm[i1]-xm[i2])/boxx));
 	   ry = ( ym[i1] - ym[i2]) - boxy* ((int)(2*(ym[i1]-ym[i2])/boxy));
 	   rz = ( zm[i1] - zm[i2]) - boxz* ((int)(2*(zm[i1]-zm[i2])/boxz));
 	   rr = sqrt(rx*rx+ry*ry+rz*rz);
	   ir = (int)(rr/grid);
         //  printf("%s\n",string_array[i1]);
	   numH10[ir]=numH10[ir]+1;
         //mon1 with mon2 and with the rest of the monomers
          for(i=start;i<end;i++){
               if(i==start+5){ 
                 //H10-2nd monomer
                       start2=start+10;
                       end2=start+10+am; 
                       for(k=start2;k<end2;k++){
   	               rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	               ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	               rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	               rr = sqrt(rx*rx+ry*ry+rz*rz);
	               ir = (int)(rr/grid);
                       //printf("%s\n",string_array[i]);
	               numH10[ir]=numH10[ir]+1;
                       }
                 //H10-all
                      start2=end2;
                      end2=(a+1)*lchain;
                      for(k=start2;k<end2;k++){
   	                rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	                ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	                rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	                rr = sqrt(rx*rx+ry*ry+rz*rz);
	                ir = (int)(rr/grid);
                        //printf("%s\n",string_array[i]);
	                 numH10[ir]=numH10[ir]+1;
                        }
                 } else if(i==start || i==(start+3) || i>(start+6)){ 
                 //O1,O4,H7,H8,H9-2nd monomer
                       start2=start+11;
                       end2=start+10+am; 
                       for(k=start2;k<end2;k++){
   	               rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	               ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	               rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	               rr = sqrt(rx*rx+ry*ry+rz*rz);
	               ir = (int)(rr/grid);
                       if(i==start) numO1[ir]=numO1[ir]+1;
                       if(i==(start+3)) numO4[ir]=numO4[ir]+1;
                       if(i>(start+6)) numH7[ir]=numH7[ir]+1; 
                       }
                 //O1,O4,H7,H8,H9-all
                      start2=end2;
                      end2=(a+1)*lchain; 
                      for(k=start2;k<end2;k++){
   	                rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	                ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	                rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	                rr = sqrt(rx*rx+ry*ry+rz*rz);
	                ir = (int)(rr/grid);
                       if(i==start) numO1[ir]=numO1[ir]+1;
                       if(i==(start+3)) numO4[ir]=numO4[ir]+1;
                       if(i>(start+6)) numH7[ir]=numH7[ir]+1; 
                       }
               } else if(i==(start+1) || i==(start+4) || i==(start+6)){
                 //C2,H5,C6-2nd monomer
                       start2=start+12;
                       end2=start+10+am; 
                       for(k=start2;k<end2;k++){
   	               rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	               ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	               rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	               rr = sqrt(rx*rx+ry*ry+rz*rz);
	               ir = (int)(rr/grid);
                       if(i==(start+1)) numC2[ir]=numC2[ir]+1;
                       if(i==(start+4)) numH5[ir]=numH5[ir]+1;
                       if(i==(start+6)) numC6[ir]=numC6[ir]+1; 
                      }
                 //C2,H5,C6-all
                      start2=end2;
                      end2=(a+1)*lchain; 
                      for(k=start2;k<end2;k++){
   	                rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	                ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	                rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	                rr = sqrt(rx*rx+ry*ry+rz*rz);
	                ir = (int)(rr/grid);
                       if(i==(start+1)) numC2[ir]=numC2[ir]+1;
                       if(i==(start+4)) numH5[ir]=numH5[ir]+1;
                       if(i==(start+6)) numC6[ir]=numC6[ir]+1; 
                      }
               } else {
                 //C3-2nd monomer
                   start2=start+14;
                   end2=start+10+am; 
                  for(k=start2;k<end2;k++){
   	             rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	             ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	             rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	             rr = sqrt(rx*rx+ry*ry+rz*rz);
	             ir = (int)(rr/grid);
	             numC3[ir]=numC3[ir]+1;
                    } 
                 //C3-all
                     start2=end2;
                     end2=(a+1)*lchain; 
                     for(k=start2;k<end2;k++){
   	              rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	              ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	              rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	              rr = sqrt(rx*rx+ry*ry+rz*rz);
	              ir = (int)(rr/grid);
	              numC3[ir]=numC3[ir]+1;
                    } 
               }
      }          
      //end of the first monomer
      //loop over inner monomers
    for(b=0;b<nmon;b++){
        start=a*lchain+b*am+10;
        end=start+am; 
    //check the exceptions within the monomer b
        i2=end-3;
          //intramonomer interactions: loop over hydrogens
        for(i=i2;i<end;i++){
         //O1
         i1=start;
   	 rx = ( xm[i1] - xm[i]) - boxx* ((int)(2*(xm[i1]-xm[i])/boxx));
 	 ry = ( ym[i1] - ym[i]) - boxy* ((int)(2*(ym[i1]-ym[i])/boxy));
 	 rz = ( zm[i1] - zm[i]) - boxz* ((int)(2*(zm[i1]-zm[i])/boxz));
 	 rr = sqrt(rx*rx+ry*ry+rz*rz);
	 ir = (int)(rr/grid);
	 numO1[ir]=numO1[ir]+1;
         //O4
         i1=start+3;
   	 rx = ( xm[i1] - xm[i]) - boxx* ((int)(2*(xm[i1]-xm[i])/boxx));
 	 ry = ( ym[i1] - ym[i]) - boxy* ((int)(2*(ym[i1]-ym[i])/boxy));
 	 rz = ( zm[i1] - zm[i]) - boxz* ((int)(2*(zm[i1]-zm[i])/boxz));
 	 rr = sqrt(rx*rx+ry*ry+rz*rz);
	 ir = (int)(rr/grid);
	 numO4[ir]=numO4[ir]+1;
         }
   //loop for mon1 mon2 interactions 
      for(i=start;i<end;i++){
          //O1,O4,H7   
         if(i==start || i==(start+3) || i>(start+5)){
                  start2=start+10;
                 if(b==(nmon-1)) end2=start+20; 
                            else end2=start+2*am;
                 for(k=start2;k<end2;k++){
   	             rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	             ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	             rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	             rr = sqrt(rx*rx+ry*ry+rz*rz);
	             ir = (int)(rr/grid);
	             if(i==start) numO1[ir]=numO1[ir]+1;
	             if(i==(start+3)) numO4[ir]=numO4[ir]+1;
	             if(i>(start+5))  numH7[ir]=numH7[ir]+1;
                 } 
                 start2=end2;
                 end2=(a+1)*lchain; 
                for(k=start2;k<end2;k++){
   	             rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	             ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	             rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	             rr = sqrt(rx*rx+ry*ry+rz*rz);
	             ir = (int)(rr/grid);
	             if(i==start) numO1[ir]=numO1[ir]+1;
	             if(i==(start+3)) numO4[ir]=numO4[ir]+1;
	             if(i>(start+5)) numH7[ir]=numH7[ir]+1;
                }
          //C2,H5,C6   
              } else if(i==(start+1) || i==(start+4) || i==(start+5)) {
                  start2=start+11;
                  if(b==(nmon-1)) end2=start+20; 
                             else end2=start+2*am;
                 for(k=start2;k<end2;k++){
   	             rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	             ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	             rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	             rr = sqrt(rx*rx+ry*ry+rz*rz);
	             ir = (int)(rr/grid);
	             if(i==(start+1))  numC2[ir]=numC2[ir]+1;
	             if(i==(start+4))  numH5[ir]=numH5[ir]+1;
	             if(i==(start+5)) numC6[ir]=numC6[ir]+1;
                 } 
                 start2=end2;
                 end2=(a+1)*lchain; 
                for(k=start2;k<end2;k++){
   	             rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	             ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	             rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	             rr = sqrt(rx*rx+ry*ry+rz*rz);
	             ir = (int)(rr/grid);
	             if(i==(start+1))  numC2[ir]=numC2[ir]+1;
	             if(i==(start+4))  numH5[ir]=numH5[ir]+1;
	             if(i==(start+5))  numC6[ir]=numC6[ir]+1;
                }
                    
               } else { 
                 //C3
                    start2=start+13;
                    if(b==(nmon-1)) end2=start+20; 
                               else end2=start+2*am;
                 for(k=start2;k<end2;k++){
   	             rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	             ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	             rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	             rr = sqrt(rx*rx+ry*ry+rz*rz);
	             ir = (int)(rr/grid);
	             numC3[ir]=numC3[ir]+1;
                  } 
                    start2=end2;
                    end2=(a+1)*lchain; 
                 for(k=start2;k<end2;k++){
   	             rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	             ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	             rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	             rr = sqrt(rx*rx+ry*ry+rz*rz);
	             ir = (int)(rr/grid);
	             numC3[ir]=numC3[ir]+1;
                 }
             }
     }
   }//end of the loop over nmon-1 monomers of the chain 
  //last monomer
    start=a*lchain+nmon*am+10;
    end=start+11; 
    //check the exceptions within the monomer
    i2=end-3;
    //if(a==nchain-1)printf("%d %d\n",start,end);
    //intramonomer interactions: loop over hydrogens
    for(i=i2;i<end;i++){
         //O1
         i1=start;
   	 rx = ( xm[i1] - xm[i]) - boxx* ((int)(2*(xm[i1]-xm[i])/boxx));
 	 ry = ( ym[i1] - ym[i]) - boxy* ((int)(2*(ym[i1]-ym[i])/boxy));
 	 rz = ( zm[i1] - zm[i]) - boxz* ((int)(2*(zm[i1]-zm[i])/boxz));
 	 rr = sqrt(rx*rx+ry*ry+rz*rz);
	 ir = (int)(rr/grid);
	 numO1[ir]=numO1[ir]+1;
         //O4
         i1=start+3;
   	 rx = ( xm[i1] - xm[i]) - boxx* ((int)(2*(xm[i1]-xm[i])/boxx));
 	 ry = ( ym[i1] - ym[i]) - boxy* ((int)(2*(ym[i1]-ym[i])/boxy));
 	 rz = ( zm[i1] - zm[i]) - boxz* ((int)(2*(zm[i1]-zm[i])/boxz));
 	 rr = sqrt(rx*rx+ry*ry+rz*rz);
	 ir = (int)(rr/grid);
	 numO4[ir]=numO4[ir]+1;
         //H11
         i1=start+6;
   	 rx = ( xm[i1] - xm[i]) - boxx* ((int)(2*(xm[i1]-xm[i])/boxx));
 	 ry = ( ym[i1] - ym[i]) - boxy* ((int)(2*(ym[i1]-ym[i])/boxy));
 	 rz = ( zm[i1] - zm[i]) - boxz* ((int)(2*(zm[i1]-zm[i])/boxz));
 	 rr = sqrt(rx*rx+ry*ry+rz*rz);
	 ir = (int)(rr/grid);
	 numH11[ir]=numH11[ir]+1;
     }
  // printf("%s %s\n",string_array[i1],string_array[i2]);
         //O1 H11
         i1=start;
         i2=start+6;
  // printf("%s %s\n",string_array[i1],string_array[i2]);
   	 rx = ( xm[i1] - xm[i2]) - boxx* ((int)(2*(xm[i1]-xm[i2])/boxx));
 	 ry = ( ym[i1] - ym[i2]) - boxy* ((int)(2*(ym[i1]-ym[i2])/boxy));
 	 rz = ( zm[i1] - zm[i2]) - boxz* ((int)(2*(zm[i1]-zm[i2])/boxz));
 	 rr = sqrt(rx*rx+ry*ry+rz*rz);
	 ir = (int)(rr/grid);
	 numO1[ir]=numO1[ir]+1;
         //O4 H11
         i1=start+3;
         i2=start+6;
       // printf("%s %s\n",string_array[i1],string_array[i2]);
   	 rx = ( xm[i1] - xm[i2]) - boxx* ((int)(2*(xm[i1]-xm[i2])/boxx));
 	 ry = ( ym[i1] - ym[i2]) - boxy* ((int)(2*(ym[i1]-ym[i2])/boxy));
 	 rz = ( zm[i1] - zm[i2]) - boxz* ((int)(2*(zm[i1]-zm[i2])/boxz));
 	 rr = sqrt(rx*rx+ry*ry+rz*rz);
	 ir = (int)(rr/grid);
	 numO4[ir]=numO4[ir]+1;
//last monomer with the rest of the atoms    
 for(i=start;i<end;i++){
         start2=end;
         end2=(a+1)*lchain; 
      //loop over the rest of the atoms 
                for(k=start2;k<end2;k++){
   	             rx = ( xm[i] - xm[k]) - boxx* ((int)(2*(xm[i]-xm[k])/boxx));
 	             ry = ( ym[i] - ym[k]) - boxy* ((int)(2*(ym[i]-ym[k])/boxy));
 	             rz = ( zm[i] - zm[k]) - boxz* ((int)(2*(zm[i]-zm[k])/boxz));
 	             rr = sqrt(rx*rx+ry*ry+rz*rz);
	             ir = (int)(rr/grid);
	             if(i==(start))  numO1[ir]=numO1[ir]+1;
	             if(i==(start+1)) numC2[ir]=numC2[ir]+1;
	             if(i==(start+2)) numC3[ir]=numC3[ir]+1;
	             if(i==(start+3)) numO4[ir]=numO4[ir]+1;
	             if(i==(start+4)) numH5[ir]=numH5[ir]+1;
	             if(i==(start+5)) numO10[ir]=numO10[ir]+1;
	             if(i==(start+6)) numH11[ir]=numH11[ir]+1;
	             if(i==(start+7)) numC6[ir]=numC6[ir]+1;
	             if(i>(start+7))  numH7[ir]=numH7[ir]+1;
                }
   }
}//end of the loop over various chains
    for(ir=0;ir<ngridwr;ir++){
        ifac=pow((ir+1),3)-pow((ir),3);
       grO1[ir]=  grO1[ir]+2*numO1[ir]/(rho*ifac*nchain*(nmon+2));
       grC2[ir]=  grC2[ir]+2*numC2[ir]/(rho*ifac*nchain*(nmon+2));
       grC3[ir]=  grC3[ir]+2*numC3[ir]/(rho*ifac*nchain*(nmon+2));
       grO4[ir]=  grO4[ir]+2*numO4[ir]/(rho*ifac*nchain*(nmon+2));
       grH5[ir]=  grH5[ir]+2*numH5[ir]/(rho*ifac*nchain*(nmon+2));
       grC6[ir]=  grC6[ir]+2*numC6[ir]/(rho*ifac*nchain*(nmon+2));
       grH7[ir]=  grH7[ir]+2*numH7[ir]/(rho*ifac*nchain*3*(nmon+2));
       grO10[ir]= grO10[ir]+2*numO10[ir]/(rho*ifac*nchain);
       grH10[ir]= grH10[ir]+2*numH10[ir]/(rho*ifac*nchain);
       grH11[ir]= grH11[ir]+2*numH11[ir]/(rho*ifac*nchain);
    }
}
 //end of the file reading
//free the memory
free(nmol);
free(natom);
			fprintf (out1,"#grO1,grC2,grC3,grO4,grH5,grC6,grH7,grO10,grH10,grH11");
	for (ir=0; ir<ngridwr; ir++){
		if (ir*grid<=boxx/2){
			fprintf (out1,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",ir*grid,grO1[ir]/(iframes),grC2[ir]/(iframes),grC3[ir]/(iframes),grO4[ir]/(iframes),grH5[ir]/(iframes),grC6[ir]/(iframes),grH7[ir]/(iframes),grO10[ir]/(iframes),grH10[ir]/(iframes),grH11[ir]/(iframes));
		}
	}
 //end of the file reading
fclose(init);
fclose(out1);
return 0;
}
