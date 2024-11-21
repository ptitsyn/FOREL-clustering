//Copyright(c) Andrey Ptitsyn, PBRC 2002
// forel.cpp : Defines the entry point for the console application.
//
#include <time.h>
#include "cluster.h"

class DistanceMatrix:public MatrixObject{
public:
	float maxdist;
	float mindist;
	float variace;
	float mean;
	float median;
	int symmetric;
	
	DistanceMatrix();
	DistanceMatrix(int);
	
	void get_stats();
	int check_symmetry();
	DistanceMatrix *Read(class Args *);
};

class Args{
public:
	int precalc;		//Is pre-calculated distance matrix available?
	int dimension;		//Number of elements in the data vector
	int metric;			//metric type
	int number;			//number of cases to cluster
   int time;
	char *datafile;		//Initial raw data file name
	char *dmatrixfile;	//Distance matrix file name, if available
	char *outputfile;	//Output file name, if not default
	char *configfile;	//Configuration file
	class DistanceMatrix *dmat;
	ClusterList *data;

	Args(int argc, char *argv[]); 

private:
	int read();
	void reminder();
};


void Args::reminder(){
	printf("\nUse:");
	printf("\nforel configfile\n");
	printf("\nConfiguration file must contain a valid set of parameters\n");
}

Args::Args(int argc, char *argv[]){
	if( argc!=2) {
		reminder();
      exit(-1);
	}else{
		configfile=argv[1];
		if(!read()) {
			reminder();
         exit(-1);
		}
	}
}

int Args::read(){
char s[100];
char *c;
FILE *f,*ff;
	f=fopen(configfile,"r");
	if(!f) return(0);
	while(fgets(s,100,f)){
		if(strstr(s,"DISTANCE_MATRIX=")){
			s[strlen(s)-1]=0;
			c=strchr(s,'=');
			ff=fopen(c+1,"r");
			if(!ff) {
				precalc=0;
			}else{
				fclose(ff);
				dmatrixfile=(char *)calloc(strlen(c),sizeof(char));
				strcpy(dmatrixfile,c+1);
				precalc=1;
			}
		}
		if(strstr(s,"DATA=")){
			s[strlen(s)-1]=0;
			c=strchr(s,'=');
			ff=fopen(c+1,"r");
			if(!ff) {
				return(0);
			}else{
				fclose(ff);
				datafile=(char *)calloc(strlen(c),sizeof(char));
				strcpy(datafile,c+1);
			}
		}
		if(strstr(s,"OUTPUT=")){
			s[strlen(s)-1]=0;
			c=strchr(s,'=');
			ff=fopen(c+1,"w");
			if(!ff) {
				outputfile="forel.out";
			}else{
				fclose(ff);
				outputfile=(char *)calloc(strlen(c),sizeof(char));
				strcpy(outputfile,c+1);
			}
		}
		if(strstr(s,"DIMENSION=")){
			c=strchr(s,'=');
			sscanf(c+1,"%d",&dimension);
		}
		if(strstr(s,"NUMBER=")){
			c=strchr(s,'=');
			sscanf(c+1,"%d",&number);
		}
		if(strstr(s,"METRIC=")){
			c=strchr(s,'=');
         c++;
			if(strstr(c,"HEMMING")) metric=1;
			if(strstr(c,"EUCLID")) metric=2;
			if(strstr(c,"CORR")) metric=3;
			if(strstr(c,"CHI2")) metric=4;
			if(strstr(c,"ZHIV")) metric=5;
			if(strstr(c,"CULB")) metric=6;
			if(strstr(c,"CITYBLOCK")) metric=7;
			if(strstr(c,"KOLM")) metric=8;
		}
	}
	return(1);
}

ClusterList *ReadText(Args *args){
int i,j;
float x;
char s[1000],*ss;
char name[1000];
ClusterPlus *cp;
ClusterElement *ce;
FILE *f;
	f=fopen(args->datafile,"r");

	if(fscanf(f,"%s",name)==EOF) return(NULL);

	if(fscanf(f,"%s",s)==EOF) return(NULL);
	ss=(char *)calloc(strlen(s),sizeof(char));
	strcpy(ss,s);
	ce= new ClusterElement(args->dimension,FLOAT);
	ce->name=ss;
	for(i=0;i<args->dimension;i++) {
   	fscanf(f,"%f",&x);
      ce->vec[i]=x;
   }
   ce->metric=args->metric;
	cp = new ClusterPlus();
	cp->Add(ce);
	ClusterList *list = new ClusterList(cp);
	j=1;
	while(j<args->number){
		if(fscanf(f,"%s",s)==EOF) {
			args->number=j;
      	break;
      }
//      if(j>4940){
//      	printf("%s\n",s);
//      }
		ss=(char *)calloc(strlen(s)+15,sizeof(char));
		strcpy(ss,s);
		ce= new ClusterElement(args->dimension,FLOAT);
		ce->name=ss;
		for(i=0;i<args->dimension;i++) {
   		fscanf(f,"%f",&x);
   	   ce->vec[i]=x;
	   }
	   ce->metric=args->metric;
		cp = new ClusterPlus();
		cp->Add(ce);
		cp->name=ss;
		list->append(cp);
		j++;
	}
	return(list);
}

void WriteTextClusters(ClusterList *list, Args *args){
int i;
ClusterList *cl;
int *d;
FILE *f;

	f=fopen(args->outputfile,"w");
	fprintf(f,"Parameters:\n");
	fprintf(f,"Configuration file:\t%s\n",args->configfile);
	fprintf(f,"Data file:\t%s\n", args->datafile);
	fprintf(f,"Dimension:\t%d\n", args->dimension);
	fprintf(f,"Number of objects:\t%d\n",args->number);
	fprintf(f,"Processing time (sec):\t%d\n",args->time);
   fprintf(f,"Cluster sizes:\n");
   d=(int *)calloc(args->number,sizeof(int));
	cl=list->go_root();
	while(1){
   	d[cl->cluster->NumberOfElements]++;
      if(!cl->next)break;
		cl=cl->go_next();
	}
   fprintf(f,"1 element(singletons):\t%d\n",d[1]);
   for(i=2;i<args->number;i++){
		if(d[i]){
		   fprintf(f,"%d elements:\t%d\n",i,d[i]);
      }
   }
   free(d);
	fprintf(f,"\n");

	fprintf(f,"\nClustering results:\n");
	cl=list->go_root();
	i=1;
	while(1){
   	cl->cluster->number=i++;
		cl->cluster->WriteText(f);
      if(!cl->next)break;
		cl=cl->go_next();
	}
	fclose(f);
}

void WriteFlatTextClusters(ClusterList *list, Args *args){
int i;
ClusterList *cl;
FILE *f;

	f=fopen(args->outputfile,"w");
	cl=list->go_root();
	i=1;
	while(1){
   	cl->cluster->number=i++;
		cl->cluster->WriteFlatText(f);
      if(!cl->next)break;
		cl=cl->go_next();
	}
	fclose(f);
}

void WriteXMLClusters(ClusterList *list, Args *args){
int i;
ClusterList *cl;
FILE *f;

	f=fopen(args->outputfile,"w");
	fprintf(f,"<PARAMETERS ");
	fprintf(f,"CONFIGFILE=%s ",args->configfile);
	fprintf(f,"DATAFILE=%s ", args->datafile);
	fprintf(f,"DIMENSION=%d ", args->dimension);
	fprintf(f,"NUMOBJECTS=%d ",args->number);
	fprintf(f,"</PARAMETERS>\n");
	fprintf(f,"<RESULTS>\n");
	cl=list->go_root();
   i=1;
	while(cl->next){
   	cl->cluster->number=i++;
		cl->cluster->WriteXML(f);
		cl=cl->go_next();
	}
	fprintf(f,"</RESULTS>\n");
	fclose(f);
};

DistanceMatrix::DistanceMatrix(){
}

DistanceMatrix::DistanceMatrix(int rows){
int i;
	numRows=numColumns=rows;
	matrix=(float **)calloc(numRows, sizeof(float *));
	for(i=0;i<numRows;i++) matrix[i]=(float *)calloc(numRows,sizeof(float));
}

DistanceMatrix *DistanceMatrix::Read(Args *args){
int i,j;
FILE *f;
	f=fopen(args->dmatrixfile,"r");
	MatrixObject *dmat = new MatrixObject(args->number);
	for(i=0;i<args->number;i++){
		for(j=0;j<args->number;j++) fscanf(f,"%f",dmat->matrix[i][j]);
	}
	fclose(f);
	return(this);
}

void DistanceMatrix::get_stats(){
int i,j,k;
VectorObject *v;
	mean=median=0;
	maxdist=mindist=matrix[0][1];
	v = new VectorObject(numRows*numColumns, FLOAT);
	k=0;
	for(i=0;i<numRows;i++){
		for(j=0;j<numColumns;j++){
			if(i==j)continue;
			if(matrix[i][j]>maxdist) maxdist=matrix[i][j];
			if(matrix[i][j]<mindist) mindist=matrix[i][j];
			mean+=matrix[i][j];
			v->vec[k++]=matrix[i][j];
		}
	}
	mean=mean/k;
	median=v->Median();
	delete v;
}

int DistanceMatrix::check_symmetry(){
int i,j;
	for(i=0;i<numRows;i++){
		for(j=i+1;j<numColumns;j++){
			if(matrix[i][j]!=matrix[j][i]) return(0);
		}
	}
	symmetric=1;
	return(symmetric);
};

void DataStats( float *maxdist, float *mindist, ClusterList *list){
int i;
ClusterPlus *blob;
ClusterElement *ce;
	blob = new ClusterPlus();
	list=list->go_root();
	while(1){
		ce=list->cluster->first;
		for(i=0;i<list->cluster->NumberOfElements;i++){
			blob->Add(ce);
			ce=ce->next;
		}
      if(!list->next) break;
		list=list->go_next();
	}
	blob->SetDistanceMetric();
	*maxdist=blob->MaxDistanceInside();
	*mindist=blob->MinDistanceInside();
	delete blob;
};

void DataScaleByZ( ClusterList *list){
int i,j;
ClusterPlus *blob;
ClusterElement *ce, *stdv;
	blob = new ClusterPlus();
	list=list->go_root();
	while(list->next){
		ce=list->cluster->first;
		for(i=0;i<list->cluster->NumberOfElements;i++){
			blob->Add(ce);
			ce=ce->next;
		}
		list=list->go_next();
	}
	blob->SetDistanceMetric();
   blob->SetCentroid();
   stdv=new ClusterElement(blob->first->number,FLOAT);
	blob->current=blob->first;
	for(i=0;i<blob->NumberOfElements;i++){
		for(j=0;j<blob->first->number;j++){
			stdv->vec[j]+=(blob->current->vec[j]-blob->centroid->vec[j])*(blob->current->vec[j]-blob->centroid->vec[j]);
		}
      blob->current=blob->current->next;
	}
	for(j=0;j<blob->first->number;j++){
		stdv->vec[j]=sqrt(stdv->vec[j]/blob->NumberOfElements);
	}

	list=list->go_root();
	while(list->next){
		ce=list->cluster->first;
		for(i=0;i<list->cluster->NumberOfElements;i++){
			for(j=0;j<stdv->number;j++){
         	ce->vec[j]=(ce->vec[j]-blob->centroid->vec[j])/stdv->vec[j];
         }
			ce=ce->next;
		}
		list=list->go_next();
	}

	delete blob;
   delete stdv;
};


MatrixObject *DataCovariation( ClusterList *list ){
int i;
ClusterPlus *blob;
ClusterElement *ce;
MatrixObject *cov;
	blob = new ClusterPlus();
	list=list->go_root();
	while(list->next){
		ce=list->cluster->first;
		for(i=0;i<list->cluster->NumberOfElements;i++){
			blob->Add(ce);
			ce=ce->next;
		}
		list=list->go_next();
	}
	blob->SetDistanceMetric();
	cov=blob->Covariation();
	delete blob;
   return(cov);
};

#define GRAVITY_FACTOR 0

ClusterElement **draft_cluster(ClusterList *start, float radius, int *step){
int i, flag, listsize;
float d;
float *m;
ClusterList *list;
ClusterPlus *clu;
ClusterElement **clip;

	listsize=0;
	list=start->go_root();
	while(list->next){
		if(!list->mark){
			list->cluster->mark=0;
   	   listsize++;
      }
		list=list->go_next();
	}
	clip=(ClusterElement **)calloc(listsize+15,sizeof(ClusterElement *));

	clu = new ClusterPlus();
   clu->NumberOfElements=0;
	clu->Add(start->cluster->first);
	clu->metric=start->cluster->metric;
	clu->SetCentroid();
   clu->SetDistanceMetric();
   clu->estimation=0;
	start->cluster->mark=1;
// allocate buffer for temporary storage of new captured points
   m=(float *)calloc(start->cluster->first->number+1,sizeof(float));
   for(i=0;i<start->cluster->first->number;i++){
   	m[i]=start->cluster->centroid->vec[i];
   }
	do{
		list=start->go_root();
      flag=0;
		while(list->next){
			if(!(list->cluster->mark+list->mark)) {
            d=clu->centroid->DistanceTo(list->cluster->centroid);
				if(d<=radius){
					clu->Add(list->cluster->first);
               for(i=0;i<start->cluster->first->number;i++){
               	m[i]+=list->cluster->first->vec[i];
               }
               list->cluster->mark=1;
            	flag++;
				}
         }
			list=list->go_next();
		}
		if(flag){
         //normal sensitivity with no regard to the gravity factor (presumed 1)
      	//clu->SetCentroid();
	      for(i=0;i<start->cluster->first->number;i++){
		     	clu->centroid->vec[i]=m[i]/(flag+1);
      	}
      	clu->estimation++;
      }
      //Gravity factor regulates the effect of previously
      //captured points on the hypersphere movement
      for(i=0;i<start->cluster->first->number;i++){
			m[i]=GRAVITY_FACTOR*clu->centroid->vec[i];
      }
	}while(flag);
   free(m);
   clu->current=clu->first;
   for(i=0;i<clu->NumberOfElements;i++){
   	clip[i]=clu->current;
      clu->current=clu->current->next;
	}
   *step=clu->estimation;
   delete clu;
	return(clip);
};

float quality(ClusterElement **clip){
int i;
float f;
static ClusterPlus cluster;
//ClusterPlus *cluster;
//	cluster=new ClusterPlus();
	cluster.metric=EUCLID;
	i=0;
   while(clip[i]){
   	cluster.Add(clip[i++]);
   }
//   f=cluster.MeanDistanceInside();
	cluster.SetCentroid();
   f=cluster.stdev();
//   f=cluster.Variance();
//	f=cluster.MeanDistanceToPoint(cluster.centroid);
  	cluster.NumberOfElements=0;
   cluster.first=0;
   delete cluster.centroid;
   cluster.centroid=NULL;
//   delete cluster;
	return(f);
};

int complete(ClusterList *list){
return(list->root->mark);
//int i;
//ClusterList *l;
//	l=list->go_root();
//	while(l->next){
//		i = l->cluster->NumberOfElements;
//		if(i==1) return(0);
//		l=l->go_next();
//	}
//	return(1);
};


ClusterList *Forelize(ClusterList *list, Args *args){
int i,steps, best_steps,loop,percent;
float maxdist,mindist,radius,best_radius,rstep,percentile_lower,percentile_upper,low,high;
float func;
float best_func;
ClusterList *buf;
ClusterPlus *cluster;
ClusterElement **cl, **best_cluster;
	percentile_lower=(float)0.2;//make it one of the parameters later
	percentile_upper=(float)0.0;//make it one of the parameters later
	if(args->precalc){
		maxdist=args->dmat->maxdist;
		mindist=args->dmat->mindist;
	}else{
		DataStats(&maxdist,&mindist,list);
	}

	low=mindist+(maxdist-mindist)*percentile_lower;
	high=mindist+(maxdist-mindist)*(1-percentile_upper);
	rstep=(maxdist-mindist)/25; //make it a parameter later
	while(1){
		list->cluster->SetCentroid(); //prepare singletons
      list->mark=0;
      if(!list->next)break;
		list=list->go_next();
	}
	loop=1;
   percent=0;
	do{
   	printf("\nLoop #%d",loop++);
		best_func=-1;
	   best_cluster=NULL;
      radius=low;
		while(radius<high){
	   	list=list->go_root();
			while(1){
				if(!list->mark) {				//skip finished cluster(s)
					cl=draft_cluster(list,radius,&steps);		//one walk of the sphere, return provisional cluster
					//func=quality(cl);
					//if(func>0) {
               //	func=1/func;				//estimate quality of the provisional cluster
               //}
               func=steps;
					if(func>best_func) {
						best_func=func;
						best_radius=radius;
                  best_steps=steps;
						if(best_cluster) free(best_cluster);
						best_cluster=cl;
					}else{
						i=0;
					   while(cl[i]){
					   	//delete(cl[i++]);
                     cl[i++]=NULL;
						}
		            free(cl);
      	      }
            }
            if(!list->next) break;
				list=list->go_next();
			}
			radius+=rstep;
		}
		//remove all members of the best cluster from the list,
		cluster = new ClusterPlus();
      cluster->NumberOfElements=0;
      cluster->metric=EUCLID;
      cluster->radius=best_radius;
      cluster->quality=best_func;
      cluster->steps=best_steps;
      cluster->name=best_cluster[0]->name;
		i=0;
		while(best_cluster[i]){
			list=list->go_root();
			while(list->next){
				if(list->cluster->first==best_cluster[i]) {
					buf=list->remove();
					list=list->go_root();
					cluster->Add(best_cluster[i]);
               delete buf;
					break;
				}else{
					list=list->go_next();
            }
			}
         i++;
		}
		//add the best cluster to the list as one cluster
		list->append(cluster);
		percent+=cluster->NumberOfElements;
      printf(" %d items classified, %f%% cumulative", cluster->NumberOfElements, (percent*100.)/args->number);
		list->last->mark=1;
      free(best_cluster);
		//special case, stop clustering if the best cluster is still a singleton
      //ATTENTION, not always aplicable
      //if(cluster->NumberOfElements<2) break;
	}while(!complete(list));
	return(list);
}

void Apprize(ClusterList *list, MatrixObject *cov, Args *args){
float d, sd;
ClusterList *cl,*cl1;
MatrixObject *icov;
FILE *f;

	f=fopen(args->outputfile,"a");
	cl=list->go_root();
	fprintf(f,"\nInter-cluster distance table (Kullback metric):\n");
	sd=0;
	while(1){
   	cl1=cl->go_next();
      while(1){
      	if((cl1->cluster->NumberOfElements>200)&&(cl->cluster->NumberOfElements>200)){
         	d=cl->cluster->KullbackTo(cl1->cluster);
	         fprintf(f,"%3.2e\t",d);
   	      sd+=d;
         }else{
         	fprintf(f,"N/A\t");
         }
         if(!cl1->next) break;
         cl1=cl1->go_next();
      }
      fprintf(f,"\n");
      if(!cl->next) break;
		cl=cl->go_next();
	}
   fprintf(f,"Metric summ:\t%3.2e",sd);

	fprintf(f,"\nInter-cluster distance table (Mahalanobis metric):\n");
	cl=list->go_root();
	cl->cluster->cov=cov;
   icov=cov->Inverse();
	sd=0;
	while(1){
   	cl1=cl->go_next();
      while(1){
      	d=cl->cluster->MahalanobisTo(cl1->cluster, icov);
         fprintf(f,"%3.2e\t",d);
         printf("%3.2e\t",d);
         sd+=d;
         if(!cl1->next) break;
         cl1=cl1->go_next();
      }
      fprintf(f,"\n");
      if(!cl->next) break;
		cl=cl->go_next();
	}
   fprintf(f,"Metric summ:\t%3.2e",sd);
   fclose(f);
}

int main(int argc, char* argv[])
{
time_t t1,t2;
DistanceMatrix *dm;
ClusterList *list;
MatrixObject *cov;
	Args *args = new Args(argc, argv);
	if(!args) return(-1);
   time(&t1);
	list = ReadText(args);
	if(args->precalc) {
		dm = new DistanceMatrix(args->number);
		dm->Read(args);
		dm->get_stats();
	}
   cov=DataCovariation(list);
	list=Forelize(list, args);
   args->time=time(&t2)-t1;
	WriteFlatTextClusters(list, args);
   printf("\nClustering done, cluster stats:\n");
	Apprize(list,cov,args);
	return 0;
}
