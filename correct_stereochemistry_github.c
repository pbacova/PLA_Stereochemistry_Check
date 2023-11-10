//This code changes the stereochemistry of the selected monomers. The input file "sequence_wrong" should contain a sequence of 1 or 0. The length of the sequence (i.e., the number of lines in the files) should correspond to the number of monomers in the system. The labeling starts from monomer 1 in the first chain (line 1) until the last monomer in the last chain (last line). 
//In "sequence_wrong" 0 is used to label monomers, which maintain their stereochemistry and 1 for monomers, whose stereochemistry will be changed by the code.
//The output is a gro file containing the same number of chains and monomers as the original file, but with the changed stereochemistry.
//Careful, the input file after_run.gro must contain unwrapped coordinates and the output file, modified.gro, contains overlaps and therefore, energy minimization procedure must be applied.
//If the input file contain also velocities, uncomment line 171 and comment the line 172. 
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include<string.h>

int nchain=70;
int lchain=70;
int mPLA=9;
void print_help()
{
	printf ("usage:re [options]\n \
	options:\n \
	-p number of chains\n \
        -k number of atoms in the chain\n \
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
                        case 'p':
                                nchain=atoi(*++argv);
                                break;
                        case 'k':
                                lchain=atoi(*++argv);
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
double r2()
{
    return (double)rand() / (double)RAND_MAX ;
}
int main(int argc, char *argv[]){
	// initialize
	double x,y,z,vx,vy,vz;
	FILE *fseq,*init,*out1;
	int  is,value1,value2,value3,value4,value5,value6,value7,value8,value9,value10,value11,value12;
        double n1,n2,n3,n4,a1,a2,a3,b1,b2,b3,t5,t6,t7,t8,t9;
        int *seq,correct,n,a,c,d,i,j,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,*nmol,*natom,start,end;
        char outname[150];
        char buf[100]; 
	char **string_array;
        double newx,newy,newz,sizea,sizeb,sizen;
        double *o1x,*o1y,*o1z,*o4x,*o4y,*o4z,*c2x,*c2y,*c2z,*c3x,*c3y,*c3z,*h5x,*h5y,*h5z,*h7x,*h7y,*h7z,*h8x,*h8y,*h8z,*h9x,*h9y,*h9z,*h10x,*h10y,*h10z,*o10x,*o10y,*o10z,*h11x,*h11y,*h11z,*c6x,*c6y,*c6z;
        double tmp7,tmp9,tmp8,decision;
        int tmp3,tmp10,tmp0; 
        char tmp1[5];
        char tempchar[5];
        parse_options (argc, argv);
	
          if((fseq=fopen("sequence_wrong","r"))==NULL){
           printf("file can not be opened\n");
           return 1;
           } 
          if((init=fopen("after_run.gro","r"))==NULL){
           printf("file can not be opened\n");
           return 1;
           } 
          int N_all=nchain*lchain;
          int mback=(lchain/mPLA);  
          sprintf(outname,"modified2.gro");
          out1=fopen(outname,"w");
          srand(time(NULL));  
          seq=(int*)malloc((mback*nchain)*sizeof(int));
           for(j=0;j<(mback*nchain);j++){
          fscanf(fseq,"%d\n",&seq[j]);
           }
          o1x=(double*)malloc((mback)*sizeof(double));
          o1y=(double*)malloc((mback)*sizeof(double));
          o1z=(double*)malloc((mback)*sizeof(double));
          o4x=(double*)malloc((mback)*sizeof(double));
          o4y=(double*)malloc((mback)*sizeof(double));
          o4z=(double*)malloc((mback)*sizeof(double));
          c2x=(double*)malloc((mback)*sizeof(double));
          c2y=(double*)malloc((mback)*sizeof(double));
          c2z=(double*)malloc((mback)*sizeof(double));
          c3x=(double*)malloc((mback)*sizeof(double));
          c3y=(double*)malloc((mback)*sizeof(double));
          c3z=(double*)malloc((mback)*sizeof(double));
          h5x=(double*)malloc((mback)*sizeof(double));
          h5y=(double*)malloc((mback)*sizeof(double));
          h5z=(double*)malloc((mback)*sizeof(double));
          h7x=(double*)malloc((mback)*sizeof(double));
          h7y=(double*)malloc((mback)*sizeof(double));
          h7z=(double*)malloc((mback)*sizeof(double));
          h8x=(double*)malloc((mback)*sizeof(double));
          h8y=(double*)malloc((mback)*sizeof(double));
          h8z=(double*)malloc((mback)*sizeof(double));
          h9x=(double*)malloc((mback)*sizeof(double));
          h9y=(double*)malloc((mback)*sizeof(double));
          h9z=(double*)malloc((mback)*sizeof(double));
          h10x=(double*)malloc((mback)*sizeof(double));
          h10y=(double*)malloc((mback)*sizeof(double));
          h10z=(double*)malloc((mback)*sizeof(double));
          o10x=(double*)malloc((mback)*sizeof(double));
          o10y=(double*)malloc((mback)*sizeof(double));
          o10z=(double*)malloc((mback)*sizeof(double));
          h11x=(double*)malloc((mback)*sizeof(double));
          h11y=(double*)malloc((mback)*sizeof(double));
          h11z=(double*)malloc((mback)*sizeof(double));
          c6x=(double*)malloc((mback)*sizeof(double));
          c6y=(double*)malloc((mback)*sizeof(double));
          c6z=(double*)malloc((mback)*sizeof(double));
          nmol=(int*)malloc(N_all*sizeof(int));
          natom=(int*)malloc(N_all*sizeof(int));
          string_array=(char**)malloc(N_all*sizeof(char*));
         for(i=0;i<N_all;i++){
          string_array[i]=(char*)malloc(5*sizeof(char));
          nmol[i]=natom[i]=0;
          }
         for(j=0;j<mback;j++)o1x[j]=o1y[j]=o1z[j]=o4x[j]=o4y[j]=o4z[j]=c2x[j]=c2y[j]=c2z[j]=c3x[j]=c3y[j]=c3z[j]=h5x[j]=h5y[j]=h5z[j]=h8x[j]=h8y[j]=h8z[j]=h9x[j]=h9y[j]=h9z[j]=c6x[j]=c6y[j]=c6z[j]=0.0; 
         for(j=0;j<mback;j++)o10x[j]=o10y[j]=o10z[j]=h11x[j]=h11y[j]=h11z[j]=0.0; 
         char str1[]="   O1"; 
         char str2[]="   C2"; 
         char str3[]="   C3"; 
         char str4[]="   O4"; 
         char str5[]="   H5"; 
         char str6[]="   C6"; 
         char str7[]="   H7"; 
         char str8[]="   H8"; 
         char str9[]="   H9"; 
         char str10[]="  H10"; 
         char str11[]="  O10"; 
         char str12[]="  H11"; 
        
        n=0;
        fgets(buf,sizeof(buf),init);
        fgets(buf,sizeof(buf),init);
           correct=0;  
        fprintf(out1,"modified file\n"); 
        fprintf(out1,"%d\n",N_all); 
    for(a=0;a<nchain;a++){
           start=a*lchain;
           end=start+lchain;
           k1=0;
           k2=0;
           k3=0;
           k4=0;
           k5=0; 
           k6=0; 
           k7=0; 
           k8=0; 
           k9=0; 
           k10=0; 
           k11=mback-1; 
           k12=mback-1; 
           for(i=start;i<end;i++){
            // fscanf(init,"%5d%5c%5c%5d%lf%lf%lf%lf%lf%lf",&nmol[n],&tmp1[0],string_array[n],&natom[n],&x,&y,&z,&vx,&vy,&vz);
             fscanf(init,"%5d%5c%5c%5d%lf%lf%lf",&nmol[n],&tmp1[0],string_array[n],&natom[n],&x,&y,&z);
             strcpy(tempchar,string_array[n]);
             value1=strcmp(tempchar,str1); 
             value2=strcmp(tempchar,str2);
             value3=strcmp(tempchar,str3);
             value4=strcmp(tempchar,str4);
             value5=strcmp(tempchar,str5); 
             value6=strcmp(tempchar,str6);
             value7=strcmp(tempchar,str7);
             value8=strcmp(tempchar,str8);
             value9=strcmp(tempchar,str9);
             value10=strcmp(tempchar,str10);
             value11=strcmp(tempchar,str11);
             value12=strcmp(tempchar,str12);
           if(value1==0){
                o1x[k1]=x;
                o1y[k1]=y;
                o1z[k1]=z;
                k1++;
                }
           if(value2==0){
                c2x[k2]=x;
                c2y[k2]=y;
                c2z[k2]=z;
                k2++;
                } 
           if(value3==0){
                c3x[k3]=x;
                c3y[k3]=y;
                c3z[k3]=z;
                k3++;
                } 
           if(value4==0){
                o4x[k4]=x;
                o4y[k4]=y;
                o4z[k4]=z;
                k4++;
                }
           if(value5==0){
                h5x[k5]=x;
                h5y[k5]=y;
                h5z[k5]=z;
                k5++;
                } 
           if(value6==0){
                c6x[k6]=x;
                c6y[k6]=y;
                c6z[k6]=z;
                k6++;
                } 
           if(value7==0){
                h7x[k7]=x;
                h7y[k7]=y;
                h7z[k7]=z;
                k7++;
                } 
           if(value8==0){
                h8x[k8]=x;
                h8y[k8]=y;
                h8z[k8]=z;
                k8++;
                } 
           if(value9==0){
                h9x[k9]=x;
                h9y[k9]=y;
                h9z[k9]=z;
                k9++;
                } 
           if(value10==0){
                h10x[k10]=x;
                h10y[k10]=y;
                h10z[k10]=z;
                k10++;
                } 
           if(value11==0){
                o10x[k11]=x;
                o10y[k11]=y;
                o10z[k11]=z;
                k11++;
                } 
           if(value12==0){
                h11x[k12]=x;
                h11y[k12]=y;
                h11z[k12]=z;
                k12++;
                } 
           n++;
           }//end of reading 1 chain
   n=n-lchain; 
  for(d=0;d<(mback-1);d++){
    a1=c2x[d]-c3x[d];
    a2=c2y[d]-c3y[d];
    a3=c2z[d]-c3z[d];
    sizea=sqrt(a1*a1+a2*a2+a3*a3);
    a1=a1/(sizea);
    a2=a2/(sizea);
    a3=a3/(sizea);
    b1=o1x[d+1]-c3x[d];
    b2=o1y[d+1]-c3y[d];
    b3=o1z[d+1]-c3z[d];
    sizeb=sqrt(b1*b1+b2*b2+b3*b3);
    b1=b1/(sizeb);
    b2=b2/(sizeb);
    b3=b3/(sizeb);
     n1=a2*b3-a3*b2;
     n2=a3*b1-a1*b3;
     n3=a1*b2-a2*b1;
     n4=-(n1*c3x[d]+n2*c3y[d]+n3*c3z[d]);
     t5=(-n4-n1*h5x[d]-n2*h5y[d]-n3*h5z[d])/(n1*n1+n2*n2+n3*n3);
     t6=(-n4-n1*c6x[d]-n2*c6y[d]-n3*c6z[d])/(n1*n1+n2*n2+n3*n3);
     t7=(-n4-n1*h7x[d]-n2*h7y[d]-n3*h7z[d])/(n1*n1+n2*n2+n3*n3);
     t8=(-n4-n1*h8x[d]-n2*h8y[d]-n3*h8z[d])/(n1*n1+n2*n2+n3*n3);
     t9=(-n4-n1*h9x[d]-n2*h9y[d]-n3*h9z[d])/(n1*n1+n2*n2+n3*n3);
     //print
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],o1x[d],o1y[d],o1z[d]); 
     n++;
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],c2x[d],c2y[d],c2z[d]); 
     n++;
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],c3x[d],c3y[d],c3z[d]); 
     n++;
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],o4x[d],o4y[d],o4z[d]); 
     n++;
//changing the chirality
     is=a*mback+d; 
  if(seq[is]==0){  
     t5=t6=t7=t8=t9=0.0;
      correct++;
    }
     newx=h5x[d]+2*t5*n1;  
     newy=h5y[d]+2*t5*n2;  
     newz=h5z[d]+2*t5*n3; 
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],newx,newy,newz); 
     n++;
     if(d==0){
       fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],h10x[d],h10y[d],h10z[d]); 
       n++;
      }
     newx=c6x[d]+2*t6*n1;  
     newy=c6y[d]+2*t6*n2;  
     newz=c6z[d]+2*t6*n3;  
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],newx,newy,newz); 
     n++;
     newx=h7x[d]+2*t7*n1;  
     newy=h7y[d]+2*t7*n2;  
     newz=h7z[d]+2*t7*n3;  
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],newx,newy,newz); 
     n++;
     newx=h8x[d]+2*t8*n1;  
     newy=h8y[d]+2*t8*n2;  
     newz=h8z[d]+2*t8*n3;  
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],newx,newy,newz); 
     n++;
     newx=h9x[d]+2*t9*n1;  
     newy=h9y[d]+2*t9*n2;  
     newz=h9z[d]+2*t9*n3;  
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],newx,newy,newz); 
     n++;
    }
     //last monom
     d=mback-1; 
     a1=c2x[d]-c3x[d];
     a2=c2y[d]-c3y[d];
     a3=c2z[d]-c3z[d];
     sizea=sqrt(a1*a1+a2*a2+a3*a3);
     a1=a1/(sizea);
     a2=a2/(sizea);
     a3=a3/(sizea);
     b1=o10x[d]-c3x[d];
     b2=o10y[d]-c3y[d];
     b3=o10z[d]-c3z[d];
     sizeb=sqrt(b1*b1+b2*b2+b3*b3);
     b1=b1/(sizeb);
     b2=b2/(sizeb);
     b3=b3/(sizeb);
     n1=a2*b3-a3*b2;
     n2=a3*b1-a1*b3;
     n3=a1*b2-a2*b1;
     n4=-(n1*c3x[d]+n2*c3y[d]+n3*c3z[d]);
     t5=(-n4-n1*h5x[d]-n2*h5y[d]-n3*h5z[d])/(n1*n1+n2*n2+n3*n3);
     t6=(-n4-n1*c6x[d]-n2*c6y[d]-n3*c6z[d])/(n1*n1+n2*n2+n3*n3);
     t7=(-n4-n1*h7x[d]-n2*h7y[d]-n3*h7z[d])/(n1*n1+n2*n2+n3*n3);
     t8=(-n4-n1*h8x[d]-n2*h8y[d]-n3*h8z[d])/(n1*n1+n2*n2+n3*n3);
     t9=(-n4-n1*h9x[d]-n2*h9y[d]-n3*h9z[d])/(n1*n1+n2*n2+n3*n3);
     //print
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],o1x[d],o1y[d],o1z[d]); 
     n++;
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],c2x[d],c2y[d],c2z[d]); 
     n++;
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],c3x[d],c3y[d],c3z[d]); 
     n++;
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],o4x[d],o4y[d],o4z[d]); 
     n++;
//changing the chirality
     is=a*mback+d; 
 if(seq[is]==0){
      t5=t6=t7=t8=t9=0.0;
      correct++;
      }  
     newx=h5x[d]+2*t5*n1;  
     newy=h5y[d]+2*t5*n2;  
     newz=h5z[d]+2*t5*n3;  
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],newx,newy,newz); 
     n++;
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],o10x[d],o10y[d],o10z[d]); 
     n++;
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],h11x[d],h11y[d],h11z[d]); 
     n++;
     newx=c6x[d]+2*t6*n1;  
     newy=c6y[d]+2*t6*n2;  
     newz=c6z[d]+2*t6*n3;  
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],newx,newy,newz); 
     n++;
     newx=h7x[d]+2*t7*n1;  
     newy=h7y[d]+2*t7*n2;  
     newz=h7z[d]+2*t7*n3;  
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],newx,newy,newz); 
     n++;
     newx=h8x[d]+2*t8*n1;  
     newy=h8y[d]+2*t8*n2;  
     newz=h8z[d]+2*t8*n3;  
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],newx,newy,newz); 
     n++;
     newx=h9x[d]+2*t9*n1;  
     newy=h9y[d]+2*t9*n2;  
     newz=h9z[d]+2*t9*n3;  
     fprintf(out1,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",a+1,"PLA",string_array[n],natom[n],newx,newy,newz); 
     n++;
  }//end of reading all chains 
        fscanf(init,"%lf%lf%lf",&tmp7,&tmp8,&tmp9);
        fprintf(out1,"%8.4f%8.4f%8.4f",tmp7,tmp8,tmp9); 
        fgets(buf,sizeof(buf),init);
   printf("number of wrong (corrected) monomers %d\n", mback*nchain-correct);    
//defining vectors 

for(i=0;i<N_all;i++)free(string_array[i]);
free(natom);
free(nmol);

fclose(init);
fclose(out1);
return 0;
}
