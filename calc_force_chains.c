/* Identify the largest cluster in each frame of a movie */
/* Revisions log:
 * 8.05.16 DM modifying to write out force chains.
 * 6.02.16 include snapshot output of particles in the system
 * 5.13.16 default input name "smtest"
 * 7.17.15 Replacing ancient Stuart/Jared lookup table with my newly
           written version.
 * 7.15.15 Find overall largest cluster also.
 * 6.19.15 Writing a read_ascii subroutine to test float->double issues
 * 9.6.13 Well, that didn't work.  Next, I try using Hermann's method.
 * 8.29.13 Written */
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

#define PI 3.14159265359

#define MAXPAR 10000
#define MAXSTACK 50000

// HARDWIRED
#define NCELLS_X 12
#define NCELLS_Y 12
#define N_CELLS 144
#define NULL_CELL N_CELLS // Empty cell representing outside the sample
#define MAXCELL 50 // Max number of grains allowed in a single cell.

struct vortex{
  int id;
  int color;
  double x;
  double y;
  double radius;
  int clusterid;
};

struct syssize{
  double SX;
  double SY;
  double SX2;
  double SY2;
};

struct mylookup{
  int num;
  int *list;
};

struct lookupdata{
  int ncellsx;
  int ncellsy;
  double xscale;
  double yscale;
};
  

void cell_interact(struct vortex *vortex,int nV,struct mylookup **lookuptable,
		   struct lookupdata lookupdata,struct syssize syssize,
		   int **distarray,
		   FILE *force_chain_out);
void distance(double *dr,double *dx,double *dy,double x1,double y1,
	      double x2,double y2,struct syssize syssize);
void get_parameters_file(double *density,struct syssize *syssize,
			 int *runtime,struct lookupdata *lookupdata,
			 int *maxcell,int *maxnum);
void map_particles_to_cells(struct vortex *vortex,int nV,
			    struct mylookup **lookuptable,
			    struct lookupdata lookupdata);
void pair_interact(int id1,int id2,struct vortex *vortex,
		   struct syssize syssize,int **distarray,
		   FILE *force_chain_out);
void parse_snapshot(int frame,int nV,struct vortex *vortex,FILE *out3);


main(int argc,char *argv[])
{
  FILE *in,*out,*out2,*out3, *force_chain_out;
  double density;
  struct syssize syssize;
  int runtime;
  int maxcell;
  struct mylookup **lookuptable;
  struct lookupdata lookupdata;
  int completeframe;
  struct vortex *vortex;
  int maxnum;
  
  char filename[120]="smtest";
  int nV,time,nVold;
  int i,j,k,id,color;
  double x,y,rad;
  int numcluster;
  double SX,SY,SX2,SY2;
  double dx,dy,dist;
  int *stack;
  int height;
  int corrected=0;
  int index,index2,marker;
  int **distarray;
  double mindist;
  int ii,jj;
  int maxcluster;
  int overallmax;
  int frame;
  int avgmax;
  int maxcount;
  int first[N_CELLS+1];
  int last[N_CELLS+1];
  int next[MAXPAR];
  int cell;
  int n,nn;
  double range2;
  int **cluster;
  int length[MAXPAR];
  int null;

  null=-1000;

  stack=(int *)malloc((int)((MAXSTACK)*sizeof(int)));
  cluster=(int **)malloc((int)((MAXPAR)*sizeof(int *)));
  cluster[0]=(int *)malloc((int)((MAXPAR*MAXPAR)*sizeof(int)));
  for(i=1;i<=MAXPAR;i++){
    cluster[i]=cluster[i-1]+MAXPAR;
  }
  distarray=(int **)malloc((int)((MAXPAR)*sizeof(int *)));
  distarray[0]=(int *)malloc((int)((MAXPAR*MAXPAR)*sizeof(int)));
  for(i=1;i<=MAXPAR;i++){
    distarray[i]=distarray[i-1]+MAXPAR;
  }

  get_parameters_file(&density,&syssize,&runtime,&lookupdata,&maxcell,&maxnum);

  vortex=malloc(maxnum*sizeof(struct vortex));
  lookuptable=(struct mylookup **)malloc(lookupdata.ncellsx
					 *sizeof(struct mylookup *));

  for(i=0;i<lookupdata.ncellsx;i++){
    lookuptable[i]=malloc(lookupdata.ncellsy*sizeof(struct mylookup));
  }
  for(j=0;j<lookupdata.ncellsy;j++){
    for(i=0;i<lookupdata.ncellsx;i++){
      lookuptable[i][j].num=0;
      lookuptable[i][j].list=malloc(maxcell*sizeof(int));
    }
  }

  if(argc<2){
    // printf("Enter name of movie: ");
    // scanf("%s",filename);
  }
  else
    sscanf(argv[1],"%s",filename);

  if((in=fopen(filename,"r"))==NULL){
    printf("Error opening file %s\n",filename);
    exit(-1);
  }

  out=fopen("bigcluster.dat","w");
  out2=fopen("lastcluster.dat","w");
  out3=fopen("snapshot.dat","w");

  //make a directory to place the files:
  mkdir("force_chain_data",0700);
  
  //create file names using string methods of c
  char str_directory[100];
  char str_time[10];

  avgmax=0;
  maxcount=0;
  frame=0;
  overallmax=0;
  while(!feof(in)){
    frame++;
    nVold=nV;
    completeframe=read_frame(in,&nV,&time,vortex);

    sprintf(str_time,"%08d",time); //convert current to a string
    strcpy(str_directory,"force_chain_data/data_t=");
    strcat(str_directory,str_time);
    force_chain_out=fopen(str_directory,"w");
    
    fprintf(force_chain_out,"#%d  %d\n", time, nV);    

    if(!completeframe) break;
    // Prepare a distance array.
    for(i=0;i<nV;i++){
      for(j=0;j<nV;j++){
	distarray[j][i]=null;
      }
    }

    map_particles_to_cells(vortex,nV,lookuptable,lookupdata);
    cell_interact(vortex,nV,lookuptable,lookupdata,syssize,distarray,force_chain_out);

    // close the current version of the file.  make a new one for a new time.
    fclose(force_chain_out);
  
    // We've now marked all particles that are touching each other.
    // Next we must assemble them into clusters.
    numcluster=0;
    // Start by making each particle into its own cluster.
    for(i=0;i<nV;i++){
      vortex[i].clusterid=i;
      cluster[i][0]=i;
      length[i]=1;
      numcluster++;
    }
    // Check the distance for all pairs of particles provided that the
    // particles are not already in the same cluster.
    for(i=0;i<nV;i++){
      for(j=0;j<nV;j++){
	if(i==j) continue;
	if(vortex[i].clusterid==vortex[j].clusterid) continue;
	// If the particles are touching, merge their clusters.
	if(distarray[i][j]>0){
	  marker=vortex[j].clusterid;
	  for(ii=0;ii<length[marker];ii++){
	    vortex[cluster[marker][ii]].clusterid=vortex[i].clusterid;
	    cluster[vortex[i].clusterid][length[vortex[i].clusterid]]=
	      cluster[marker][ii];
	    (length[vortex[i].clusterid])++;
	  }
	  length[marker]=0;
	  numcluster--;
	}
      }
    }
	

    //printf("found %d clusters\n",numcluster);
      
    // Debug: write out results.
    maxcluster=0;
    for(i=0;i<nV;i++){
      for(j=0;j<length[i];j++){
	if(length[i]>maxcluster) maxcluster=length[i];
      }
    }
    fprintf(out,"%d %d\n",frame,maxcluster);
    avgmax+=maxcluster;
    maxcount++;
    if(overallmax<maxcluster) overallmax=maxcluster;
    if(!(frame%500))
    {
      parse_snapshot(frame,nV,vortex,out3);
      printf("Frame %d\n",frame);
    }
  }
  printf("Average maximum cluster size: %e\n",(double)avgmax/(double)maxcount);
  printf("Scaled maximum cluster size: %e\n",(double)avgmax/(double)(maxcount*nV));
  printf("Absolute maximum cluster size: %d\n",overallmax);
  printf("Reached frame %d\n",frame);
  // Print out final frame
  for(i=0;i<nV;i++){
    for(j=0;j<length[i];j++){
      fprintf(out2,"%e %e %d\n",vortex[cluster[i][j]].x,
	      vortex[cluster[i][j]].y,i);
    }
    fprintf(out2,"\n");
  }
  fclose(out);
  fclose(out2);
  fclose(out3);

}
int read_frame(FILE *in,int *nV,int *time,struct vortex *vortex)
{
  int i;
  float xin,yin,zin;

  fread(nV,sizeof(int),1,in);
  fread(time,sizeof(int),1,in);
  if(feof(in)) return 0;
  for(i=0;i<*nV;i++){
    fread(&(vortex[i].color),sizeof(int),1,in);
    fread(&(vortex[i].id),sizeof(int),1,in);
    fread(&xin,sizeof(float),1,in);
    fread(&yin,sizeof(float),1,in);
    fread(&zin,sizeof(float),1,in);
    vortex[i].x=(double)xin;
    vortex[i].y=(double)yin;
    vortex[i].radius=(double)zin;
    vortex[i].clusterid=-1;
  }
  return 1;
  
}
void cell_interact(struct vortex *vortex,int nV,struct mylookup **lookuptable,
		   struct lookupdata lookupdata,struct syssize syssize,
		   int **distarray, FILE *force_chain_out)
{
  int i,j,k,kk;
  int id1,id2;
  int loop;
  int indx,jndx;

  for(i=0;i<lookupdata.ncellsx;i++){
    for(j=0;j<lookupdata.ncellsy;j++){
      for(k=0;k<lookuptable[i][j].num;k++){
	id1=lookuptable[i][j].list[k];
	for(loop=0;loop<5;loop++){
	  switch(loop){
	  case 0:
	    // Same cell
	    indx=i;
	    jndx=j;
	    break;
	  case 1:
	    // Cell to the right
	    indx=i+1;
	    jndx=j;
	    if(indx>=lookupdata.ncellsx) indx-=lookupdata.ncellsx;
	    break;
	  case 2:
	    // Cell to diagonal upper right
	    indx=i+1;
	    jndx=j+1;
	    if(indx>=lookupdata.ncellsx) indx-=lookupdata.ncellsx;
	    if(jndx>=lookupdata.ncellsy) jndx-=lookupdata.ncellsy;
	    break;
	  case 3:
	    // Cell above
	    indx=i;
	    jndx=j+1;
	    if(jndx>=lookupdata.ncellsy) jndx-=lookupdata.ncellsy;
	    break;
	  case 4:
	    // Cell to diagonal upper left
	    indx=i-1;
	    jndx=j+1;
	    if(indx<0) indx+=lookupdata.ncellsx;
	    if(jndx>=lookupdata.ncellsy) jndx-=lookupdata.ncellsy;
	    break;
	  }
	  for(kk=0;kk<lookuptable[indx][jndx].num;kk++){
	    // Do each interacting pair only once
	    if((!loop)&&(kk<=k)) continue;
	    id2=lookuptable[indx][jndx].list[kk];
	    pair_interact(id1,id2,vortex,syssize,distarray,force_chain_out);
	  }
	}
      }
    }
  }
}
void pair_interact(int id1,int id2,struct vortex *vortex,
		   struct syssize syssize,int **distarray,
		   FILE *force_chain_out)
{
  double dr,dx,dy;
  double deltar;
  double epsilon;

  epsilon=5e-2;
  
  distance(&dr,&dx,&dy,vortex[id1].x,vortex[id1].y,vortex[id2].x,
	   vortex[id2].y,syssize);
  deltar=dr-(vortex[id1].radius+vortex[id2].radius);
  if(deltar<epsilon){
    distarray[id1][id2]=1;
    distarray[id2][id1]=1;

    //DM print the overlap vector and particle information
    //id1 x1 y1 radius1 color1 
    fprintf(force_chain_out,"%d  %f  %f  %f  %d  ", id1, vortex[id1].x,vortex[id1].y, vortex[id1].radius, vortex[id1].color);
    //id2 x2 y2 radius2 color2 
    fprintf(force_chain_out,"%d  %f  %f  %f  %d  ", id2, vortex[id2].x,vortex[id2].y, vortex[id2].radius, vortex[id2].color);
    //deltar dr dx dy
    fprintf(force_chain_out,"%f  %f  %f  %f\n", deltar, dr, dx, dy);    
  }
}
void map_particles_to_cells(struct vortex *vortex,int nV,
			    struct mylookup **lookuptable,
			    struct lookupdata lookupdata)
{
  int i,j,k;
  int xindex,yindex;

  for(j=0;j<lookupdata.ncellsy;j++){
    for(i=0;i<lookupdata.ncellsx;i++){
      lookuptable[i][j].num=0;
    }
  }

  for(i=0;i<nV;i++){
    xindex=floor(vortex[i].x*lookupdata.xscale);
    yindex=floor(vortex[i].y*lookupdata.yscale);
    if(xindex==lookupdata.ncellsx)
      xindex=0;
    if(yindex==lookupdata.ncellsy)
      yindex=0;
    k=lookuptable[xindex][yindex].num;
    lookuptable[xindex][yindex].list[k]=i;
    (lookuptable[xindex][yindex].num)++;
  }
}
void get_parameters_file(double *density,struct syssize *syssize,
			 int *runtime,struct lookupdata *lookupdata,
			 int *maxcell,int *maxnum)
{
  FILE *in;
  char trash[120];
  double pdensity;
  double cellsize;
  double resolution;
  double radius;
  double runforce;
  double dt;
  int maxtime;
  int writemovietime;
  double kspring;

  resolution=1e-6;

  if((in=fopen("Pa0","r"))==NULL){
    printf("Input file Pa0 not found\n");
    exit(-1);
  }
  fscanf(in,"%s %lf\n",trash,density);
  fscanf(in,"%s %lf\n",trash,&pdensity);
  fscanf(in,"%s %lf\n",trash,&((*syssize).SX));
  (*syssize).SX2=(*syssize).SX*0.5;
  fscanf(in,"%s %lf\n",trash,&((*syssize).SY));
  (*syssize).SY2=(*syssize).SY*0.5;
  *maxnum=*density*(*syssize).SX*(*syssize).SY;
  fscanf(in,"%s %lf\n",trash,&radius);
  fscanf(in,"%s %d\n",trash,runtime);
  fscanf(in,"%s %lf\n",trash,&runforce);
  fscanf(in,"%s %lf\n",trash,&dt);
  fscanf(in,"%s %d\n",trash,&maxtime);
  fscanf(in,"%s %d\n",trash,&writemovietime);
  fscanf(in,"%s %lf\n",trash,&kspring);
  fscanf(in,"%s %lf\n",trash,&cellsize); // size of lookup cell
  (*lookupdata).ncellsx=(int)(*syssize).SX/cellsize;
  (*lookupdata).ncellsy=(int)(*syssize).SY/cellsize;
  if(fabs((*lookupdata).ncellsx*cellsize-(*syssize).SX)>resolution){
    printf("Error, SX mismatch with cell size\n");
    exit(-1);
  }
  if(fabs((*lookupdata).ncellsy*cellsize-(*syssize).SY)>resolution){
    printf("Error, SY mismatch with cell size\n");
    exit(-1);
  }
  (*lookupdata).xscale=(double)(*lookupdata).ncellsx/(*syssize).SX;
  (*lookupdata).yscale=(double)(*lookupdata).ncellsy/(*syssize).SY;
  // Compute largest number of particles that could possibly fit
  // inside a cell of this size
  *maxcell=2*(int)cellsize*cellsize/(PI*radius*radius);
  fclose(in);
}
void distance(double *dr,double *dx,double *dy,double x1,double y1,
	      double x2,double y2,struct syssize syssize)
{
  double locdx,locdy;

  locdx=x1-x2;
  if(locdx>syssize.SX2) locdx-=syssize.SX;
  if(locdx<=-syssize.SX2) locdx+=syssize.SX;
  locdy=y1-y2;
  if(locdy>syssize.SY2) locdy-=syssize.SY;
  if(locdy<=-syssize.SY2) locdy+=syssize.SY;
  *dr=sqrt(locdx*locdx+locdy*locdy);
  *dx=locdx;
  *dy=locdy;
}
void parse_snapshot(int frame,int nV,struct vortex *vortex,FILE *out3)
{
  int i;
  fprintf(out3,"\"Frame=%d\"\n",frame);
  for(i=0;i<nV;i++)
    fprintf(out3,"%lf %lf %lf %d\n",vortex[i].x,vortex[i].y,vortex[i].radius,vortex[i].clusterid);
  fprintf(out3,"\n\n\n");
  fflush(out3);
}
