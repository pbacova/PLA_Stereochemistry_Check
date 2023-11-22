//This code calculates the internal distances of PLA chains. Check 100LA.itp file in stereochemistry_check directory for the order of the atoms in the gro file. 
////Careful, the input file must contain unwrapped coordinates.
////If the input file contain also velocities, uncomment line 122 and comment the line 123.
////The output is a dependence of n (index difference between two backbone atoms) vs. Rn^2/n, where Rn is the distance between the given two atoms i and i+n on the backbone.
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include<string.h>

int nchain=70;
int nframes=76;
int lchain=70;
int mPLA=9;
int skip=0;
char **p;
void print_help()
{
	printf ("usage:re [options]\n \
	options:\n \
	-f number of frames\n \
	-p number of chains\n \
        -k number of monomers in the chain\n \
        -o number of frames to skip\n \
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
	double x,y,z;
	FILE *init,*out2;
	int  d,n,b,i,j,k,a,c,*nmol,*natom,s,start,end;
        char outname[150];
        char buf[100]; 
	char **string_array;
        long int *nPLA;
        double *rePLA,vx,vy,vz,sumx,sumy,sumz,*xm,*ym,*zm;
        double rx,ry,rz,tmp7,tmp9,tmp8;
        int tmp3,tmp10,tmp0; 
        char tmp1[5];
        parse_options (argc, argv);
	
          if((init=fopen("input.gro","r"))==NULL){
           printf("file can not be opened\n");
           return 1;
           } 
          int iframes=nframes-skip; 
          int N_all=nchain*lchain;
          //number of backbone atoms
          int mback=3*(lchain/mPLA);  
          sprintf(outname,"in_d_all");
          out2=fopen(outname,"w");
          rePLA=(double*)malloc((mback)*sizeof(double));
          nPLA=(long int*)malloc((mback)*sizeof(long int));
          xm=(double*)malloc((mback)*sizeof(double));
          ym=(double*)malloc((mback)*sizeof(double));
          zm=(double*)malloc((mback)*sizeof(double));
          nmol=(int*)malloc(N_all*sizeof(int));
          natom=(int*)malloc(N_all*sizeof(int));
          string_array=(char**)malloc(N_all*sizeof(char*));
         for(i=0;i<N_all;i++){
         string_array[i]=(char*)malloc(6*sizeof(char));
          nmol[i]=natom[i]=0;
          }
         for(j=0;j<mback;j++)rePLA[j]=0.0; 
         for(j=0;j<mback;j++)nPLA[j]=0; 
         for(j=0;j<mback;j++)xm[j]=ym[j]=zm[j]=0.0; 
for(j=0;j<skip;j++){
        fgets(buf,sizeof(buf),init);
        fgets(buf,sizeof(buf),init);
        for(i=0;i<N_all;i++){
         fgets(buf,sizeof(buf),init);
         }
        fscanf(init,"%lf%lf%lf",&tmp7,&tmp8,&tmp9);
        fgets(buf,sizeof(buf),init);
}
for(j=0;j<(iframes);j++){
        n=0;
        fgets(buf,sizeof(buf),init);
        fgets(buf,sizeof(buf),init);
     for(a=0;a<nchain;a++){ 
           sumx=sumy=sumz=0.0;
           start=a*lchain;
           end=start+lchain;
           k=0;
           for(i=start;i<end;i++){
              fscanf(init,"%5d%5c%5c%5d%lf%lf%lf",&nmol[n],&tmp1[0],string_array[n],&natom[n],&x,&y,&z);
             // fscanf(init,"%5d%5c%5c%5d%lf%lf%lf%lf%lf%lf",&nmol[n],&tmp1[0],string_array[n],&natom[n],&x,&y,&z,&vx,&vy,&vz);
              if((i-start)<10){
                if((i-start)==0){
                xm[k]=x;
                ym[k]=y;
                zm[k]=z;
                k++;
                }
                if((i-start)==1){
                xm[k]=x;
                ym[k]=y;
                zm[k]=z;
                k++;
                }
                if((i-start)==2){
                xm[k]=x;
                ym[k]=y;
                zm[k]=z;
                k++;
                }
              }else if((i-start)<(lchain-12)){     
                if((i-start)%mPLA==1){
                xm[k]=x;
                ym[k]=y;
                zm[k]=z;
                k++;
                }
                if((i-start)%mPLA==2){
                xm[k]=x;
                ym[k]=y;
                zm[k]=z;
                k++;
                }
                if((i-start)%mPLA==3){
                xm[k]=x;
                ym[k]=y;
                zm[k]=z;
                k++;
                }
              } else {
                if((i-start)==(lchain-11)){
                xm[k]=x;
                ym[k]=y;
                zm[k]=z;
                k++;
                }
                if((i-start)==(lchain-10)){
                xm[k]=x;
                ym[k]=y;
                zm[k]=z;
                k++;
                }
                if((i-start)==(lchain-9)){
                xm[k]=x;
                ym[k]=y;
                zm[k]=z;
                k++;
                }
              }
              n++;
            }
       for(s=1;s<(mback-1);s++){
            for(d=0;d<(mback-s);d++){
               rx=xm[d+s]-xm[d]; 
               ry=ym[d+s]-ym[d]; 
               rz=zm[d+s]-zm[d];
               rePLA[s]=rePLA[s]+(rx*rx+ry*ry+rz*rz); 
               nPLA[s]++; 
               }
           }

    }//end of chains 
        fscanf(init,"%lf%lf%lf",&tmp8,&tmp8,&tmp9);
        fgets(buf,sizeof(buf),init);
   }
 //free the memory
      for(s=1;s<(mback-1);s++){
     fprintf(out2,"%d %lf\n",s,rePLA[s]/(nPLA[s]*s));
       }
fclose(init);
fclose(out2);
return 0;
}
