//Copyright(c) Andrey Ptitsyn, PBRC 2002
#include "cluster.h"

VectorObject::VectorObject(){
}

VectorObject::VectorObject(int num, int _vartype){
	number=num;                                       
	var_type=_vartype;
   name=NULL; vec=NULL; ivec=NULL;
	if(var_type==FLOAT){
		vec=(float *)calloc(number,sizeof(float));
	}
	if(var_type==INT){
		ivec=(int *)calloc(number,sizeof(int));
	}
}

VectorObject::~VectorObject(){
	if(name) {
   	free(name);
      name=NULL;
   }
	if(vec) {
   	free(vec);
      vec=NULL;
   }
	if(ivec) {
   	free(ivec);
      ivec=NULL;
   }
}

void VectorObject::Read(FILE *f){
int i;
	for(i=0;i<number;i++){
		if(var_type==FLOAT){
			fscanf(f,"%f",&vec[i]);
		}else{
			fscanf(f,"%d",&ivec[i]);
		}
	}
}

void VectorObject::Write(FILE *f){
int i;
	for(i=0;i<number;i++){
		if(var_type==FLOAT){
			fprintf(f,"%f\t",vec[i]);
		}else{
			fprintf(f,"%d\t",ivec[i]);
		}
	}
	fprintf(f,"\n");
}

float VectorObject::Mean(){
int i;
	mean=0;
	if(var_type==FLOAT){
		for(i=0;i<number;i++){
			mean+=vec[i];
		}
	}else{
		for(i=0;i<number;i++){
			mean+=ivec[i];
		}
	}
	mean/=number;
	return(mean);
}

float VectorObject::Max(){
int i;
float max;
	if(var_type==FLOAT){
		max=vec[0];
		for(i=0;i<number;i++){
			if(max<vec[i]) max=vec[i];
		}
	}else{
		max=ivec[0];
		for(i=0;i<number;i++){
			if(max<ivec[i]) max=ivec[i];
		}
	}
	return(max);
}

float VectorObject::Min(){
int i;
float min;
	if(var_type==FLOAT){
		min=vec[0];
		for(i=0;i<number;i++){
			if(min>vec[i]) min=vec[i];
		}
	}else{
   	min=ivec[0];
		for(i=0;i<number;i++){
			if(min>ivec[i]) min=ivec[i];
		}
	}
	return(min);
}

static int ISorter(const void *x, const void *y){
int *a,*b;
	a=(int *)x; b=(int *)y;
	return(*a-*b);
}

static int FSorter( const void *x, const void *y)
{
float *a,*b;
	a=(float *)x; b=(float *)y;
	return((int)(*a-*b));
}

float VectorObject::Median(){
int i;
float *buf;
int *ibuf;
	if(var_type==FLOAT){
		buf=(float *)malloc(number*sizeof(float));
		memcpy(buf,vec,number*sizeof(float));
		qsort((void *)buf,number, sizeof(float),FSorter);
		i=number%2;
		if(i){
			median=buf[i];
		}else{
			median=(buf[i]+buf[i+1])/2;
		}
		free(buf);
	}else{
		ibuf=(int *)malloc(number*sizeof(int));
		memcpy(ibuf,ivec,number*sizeof(int));
		qsort(ibuf,number, sizeof(int),ISorter);
		i=number%2;
		if(i){
			median=(float)ibuf[i];
		}else{
			median=(float)((ibuf[i]+ibuf[i+1])/2);
		}
		free(ibuf);
	}
	return(median);
}

/////////////////////////////////////////////////////////////////////

ClusterElement::ClusterElement(int num, int _vartype){
	number=num;
	var_type=_vartype;
	if(var_type==FLOAT){
		vec=(float *)calloc(number,sizeof(float));
      ivec=NULL;
	}
	if(var_type==INT){
		ivec=(int *)calloc(number,sizeof(int));
      vec=NULL;
	}
}

ClusterElement::ClusterElement(ClusterElement *ce){
int i;
	number=ce->number;
	var_type=ce->var_type;
	if(var_type==FLOAT){
		vec=(float *)calloc(ce->number,sizeof(float));
	   for(i=0;i<number;i++) vec[i]=ce->vec[i];
      ivec=NULL;
	}
	if(var_type==INT){
		ivec=(int *)calloc(ce->number,sizeof(int));
	   for(i=0;i<number;i++) ivec[i]=ce->ivec[i];
      vec=NULL;
	}
}

int ClusterElement::BinaryDistanceTo( ClusterElement *element ){
int i;
int summ=0;
if(var_type==FLOAT){
	for(i=0;i<number;i++){
		if(vec[i]!=element->vec[i]) summ++;
	}
}else{
	for(i=0;i<number;i++){
		if(ivec[i]!=element->ivec[i]) summ++;
	}
}
return(summ);
};

float ClusterElement::EuclideanDistanceTo( ClusterElement *element ){
float summ=0;
int i;
	for(i=0;i<number;i++){
		summ+=(vec[i]-element->vec[i])*(vec[i]-element->vec[i]);
	}
	summ=(float)sqrt(summ);
	return(summ);
};

float ClusterElement::CorrelationTo( ClusterElement *element ){
int i;
float m1,m2, a,b,x1,x2,x3;
	m1=this->Mean();
   m2=element->Mean();
   x1=x2=x3=0;
	for(i=0;i<number;i++){
   	a=(vec[i]-m1);
      b=(element->vec[i]-m2);
      x1+=a*b;
      x2+=a*a;
      x3+=b*b;
   }
   x2=x2*x3;
   if(x2){
	   return(1-(x1/sqrt(x2))*(x1/sqrt(x2)));
   }else{
   	return(1);
   }
}

int ClusterElement::EqualsTo( ClusterElement *element ){
int i;
	for(i=0;i<number;i++){
		if(vec[i]!=element->vec[i]) return(0);
	}
	return(1);
};

void ClusterElement::SetDistanceMetric(int DistanceMetric)
{
	metric=DistanceMetric;
}

float ClusterElement::DistanceTo(ClusterElement *another)
{
	switch(metric){
	case HEMMING:
		return( (float)BinaryDistanceTo(another));
	case EUCLID:
		return(EuclideanDistanceTo(another));
	case CORR:
		return(CorrelationTo(another));
	default: return(-1);
	}
}

////////////////////////////////////////////////////////////////////////////////

Cluster::Cluster(){
	metric=EUCLID;
   name=NULL;
   number=0;
   NumberOfElements=0;
   redundant=0;
   current=first=last=NULL;
}


void Cluster::Add(ClusterElement *element){
	if(NumberOfElements!=0){
		last->next=element;
		element->previous=last;
		element->next=NULL;
		last=element;
	}else{
		first=last=current=element;
		first->previous=NULL;
		first->next=NULL;
	}
	NumberOfElements++;
};

void Cluster::Remove(ClusterElement *element){
	if(NumberOfElements > 1){
		element->previous->next=element->next;
		element->next->previous=element->previous;
	}else{
		first=NULL;
	}
	NumberOfElements--;
};

Cluster *Cluster::Divide(){
int i;
Cluster *subcluster;
	subcluster = new Cluster;
	current=first;
	for(i=0;i<NumberOfElements;i++){
		if(current->mark!=0){
			Remove(current);
			subcluster->Add(current);
		}
		current=current->next;
	}
	return(subcluster);
};

void Cluster::Glue(Cluster *cluster){
int i;
	current=cluster->first;
	for(i=0;i<cluster->NumberOfElements;i++){
		//cluster->Remove(current);
		Add(current);
		current=current->next;
	}
};

void Cluster::Add(Cluster *cluster){
ClusterElement *ce;
int i;
	current=cluster->first;
	for(i=0;i<cluster->NumberOfElements;i++){
		ce=new ClusterElement(current);
		Add(ce);
		current=current->next;
	}
};


void Cluster::PickElementsFrom(Cluster *cluster){
int i;
	current=cluster->first;
	for(i=0;i<cluster->NumberOfElements;i++){
		if(current->mark!=0){
			Remove(current);
			Add(current);
		}
		current=current->next;
	}
};

void Cluster::GetElement(ClusterElement *element){
	Remove(element);
	Add(element);
};

int Cluster::Redundant(){
int i;
ClusterElement *element;
	current=first;
	for(i=0;i<NumberOfElements;i++){
		element=current;
		while(element->next!=NULL){
			if(current->EqualsTo(element)) return(1);
			element=element->next;
		}
		current=current->next;
	}
	return(0);
};

void Cluster::MarkAll(){
int i;
	current=first;
	for(i=0;i<NumberOfElements;i++){
		current->mark++;
		current=current->next;
	}
}


void Cluster::UnmarkAll(){
int i;
	current=first;
	for(i=0;i<NumberOfElements;i++){
		current->mark=0;
		current=current->next;
	}
}

void Cluster::SetDistanceMetric(){
	metric=first->metric;
}



/////////////////////////////////////////////////////////////////////////////
ClusterPlus::ClusterPlus(){
	estimation=0;
   steps=0;
   mark=0;
   quality=0;
   var=0;
   radius=0;
   centroid=NULL;
   cov=NULL;
   cor=NULL;
}

ClusterPlus::~ClusterPlus(){
	if(centroid) delete centroid;
   if(cov) delete cov;
   if(cor) delete cor;
}

float ClusterPlus::DistanceInside()
{
	if(estimation==MEDIAN) return(MedianDistanceInside());
	if(estimation==MEAN) return(MeanDistanceInside());
	return(0);
}

float ClusterPlus::DistanceBetween(ClusterPlus *another)
{
	if(estimation==MEDIAN) return(MedianDistanceBetween(another));
	if(estimation==MEAN) return(MeanDistanceBetween(another));
	return(0);
}

float ClusterPlus::DistanceToPoint(ClusterElement *point)
{
	if(estimation==MEDIAN) return(MedianDistanceToPoint(point));
	if(estimation==MEAN) return(MeanDistanceToPoint(point));
	return(0);
}


ClusterElement *ClusterPlus::MeanCentroid(){
int i,j;
	centroid = new ClusterElement(first->number,FLOAT);
	current=first;
	for(i=0;i<NumberOfElements;i++){
		for(j=0;j<current->number;j++){
			centroid->vec[j]+=current->vec[j];
		}
	}
	j=centroid->number/2;
	for(i=0;i<NumberOfElements;i++){
		if(centroid->vec[i]>=j){
			centroid->vec[i]=1;
		}else{
			centroid->vec[i]=0;
		}
	}
	return(centroid);
}

ClusterElement *ClusterPlus::MedianCentroid(){
int i,j;
	centroid = new ClusterElement(first->number,FLOAT);
	centroid->vec=(float *)calloc(first->number, sizeof(int)),
	current=first;
	for(i=0;i<NumberOfElements;i++){
		for(j=0;j<current->number;j++){
			centroid->vec[j]+=current->vec[j];
		}
	}
	j=centroid->number/2;
	for(i=0;i<NumberOfElements;i++){
		if(centroid->vec[i]>=j){
			centroid->vec[i]=1;
		}else{
			centroid->vec[i]=0;
		}
	}
	return(centroid);
}

ClusterElement *ClusterPlus::SetCentroid(){
int i,j;
	if(!centroid) centroid = new ClusterElement(first->number,FLOAT);
	current=first;
	for(i=0;i<NumberOfElements;i++){
		for(j=0;j<first->number;j++){
			centroid->vec[j]+=current->vec[j];
         current->mark=0;

		}
      current=current->next;
	}
	for(j=0;j<first->number;j++){
		centroid->vec[j]=centroid->vec[j]/NumberOfElements;
	}
	centroid->number=first->number;
   centroid->var_type=first->var_type;
   centroid->SetDistanceMetric(first->metric);
   centroid->name=NULL;
	return(centroid);
}

float ClusterPlus::Radius(){
int i;
float m;
VectorObject *v;
	if(!this->centroid) this->SetCentroid();
	v = new VectorObject(NumberOfElements,FLOAT);
	current=first;
	for(i=0;i<NumberOfElements;i++){
		v->vec[i]=current->DistanceTo(this->centroid);
		current=current->next;
	}
	this->radius=v->Max();
	delete v;
	return(this->radius);
}

float ClusterPlus::MinDistanceInside(){
int i,j;
float min,x;
ClusterElement *c;
	current=first; min=MAXFLOAT;
	for(i=0;i<NumberOfElements;i++){
		c=first;
		for(j=0;j<NumberOfElements;j++){
			if(current!=c) {
				x=current->DistanceTo(c);
				if(min>x) min=x;
         }
         c=c->next;
		}
		current=current->next;
	}
	return(min);
}

float ClusterPlus::MaxDistanceInside(){
int i,j;
float max,x;
ClusterElement *c;
	current=first; max=0;
	for(i=0;i<NumberOfElements;i++){
		c=first;
		for(j=0;j<NumberOfElements;j++){
			if(current!=c){
				x=current->DistanceTo(c);
				if(max<x) max=x;
         }
         c=c->next;
		}
		current=current->next;
	}
	return(max);
}

float ClusterPlus::MeanDistanceInside(){
int i,j;
float mean;
ClusterElement *c;
	current=first; mean=0;
	for(i=0;i<NumberOfElements;i++){
		c=first;
		for(j=0;j<NumberOfElements;j++){
			if(i==j) continue;
			mean+=current->DistanceTo(c);
         c=c->next;
		}
		current=current->next;
	}
	mean/=NumberOfElements;
	return(mean);
}

float ClusterPlus::MeanDistanceToPoint( ClusterElement *element){
int i;
float mean;
	mean=0;
	current=first;
	for(i=0;i<NumberOfElements;i++){
		mean+=current->DistanceTo(element);
		current=current->next;
	}
	mean/=NumberOfElements;
	return(mean);
}

float ClusterPlus::MeanDistanceBetween( ClusterPlus *another){
int i,j;
float mean;
	mean=0;
	current=first;
	for(i=0;i<NumberOfElements;i++){
		another->current=another->first;
		for(j=0;j<another->NumberOfElements;j++){
			mean+=current->DistanceTo( another->current );
			another->current=another->current->next;
		}
		current->next;
	}
	mean=mean/(NumberOfElements+another->NumberOfElements);
	return(mean);
}

float ClusterPlus::Variance(){
int i;
float var=0;
	if(NumberOfElements==1) return(var);
   current=first;
	for(i=0;i<NumberOfElements;i++){
		var+=centroid->DistanceTo(current);
		current->next;
	}
	var=var/(NumberOfElements-1);
	return(var);
}

float ClusterPlus::stdev(){
int i;
float var=0;
	if(NumberOfElements==1) return(var);
   current=first;
	for(i=0;i<NumberOfElements;i++){
		var+=centroid->DistanceTo(current);
		current->next;
	}
	var=(float)sqrt(var/(NumberOfElements-1));
	return(var);
}

float ClusterPlus::MedianDistanceInside(){
int i,j,k;
float m;
ClusterElement *c;
VectorObject *v;
	i=(NumberOfElements*NumberOfElements-NumberOfElements)/2;
	v = new VectorObject(i, FLOAT);
	current=first; k=0;
	for(i=0;i<NumberOfElements;i++){
		c=first;
		for(j=i;j<NumberOfElements;j++){
			v->vec[k++]=current->DistanceTo(c);
			c=c->next;
		}
		current=current->next;
	}
	m=v->Median();
	delete v;
	return(m);
}

float ClusterPlus::MedianDistanceToPoint( ClusterElement *element){
int i;
float m;
VectorObject *v;
	v = new VectorObject(NumberOfElements,FLOAT);
	current=first;
	for(i=0;i<NumberOfElements;i++){
		v->vec[i]=current->DistanceTo(element);
		current=current->next;
	}
	m=v->Median();
	delete v;
	return(m);
}

float ClusterPlus::MedianDistanceBetween( ClusterPlus *another){
int i,j,k;
float m;
VectorObject *v;
	v = new VectorObject((NumberOfElements*another->NumberOfElements-NumberOfElements)/2, FLOAT);
	current=first; k=0;
	for(i=0;i<NumberOfElements;i++){
		another->current=another->first;
		for(j=i;j<another->NumberOfElements;j++){
			v->vec[k++]=current->DistanceTo(another->current);
			another->current=another->current->next;
		}
		current=current->next;
	}
	m=v->Median();
	delete v;
	return(m);
}


MatrixObject *ClusterPlus::Covariation(){
	if(estimation==MEDIAN) return(MedianCovariation());
	if(estimation==MEAN) return(MeanCovariation());
	return(NULL);
}


MatrixObject *ClusterPlus::MedianCovariation(){
int i,j,k,l;
float *m;
VectorObject *buf;
MatrixObject *covariation;
	l=first->number;
	covariation= new MatrixObject( l );
	m=(float *)calloc(l,sizeof(float));
	buf= new VectorObject(NumberOfElements, FLOAT);
	for(i=0;i<l;i++){
		current=first;
		for(j=0;j<NumberOfElements;j++){
			buf->vec[j]=current->vec[i];
			current=current->next;
		}
		m[i]=buf->Median();
	}
	for(i=0;i<l;i++){
		for(j=0;j<l;j++){
			current=first;
			for(k=0;k<NumberOfElements;k++){
				buf->vec[k]=(current->vec[i]-m[i])*(current->vec[j]-m[j]);
				current=current->next;
			}
			covariation->matrix[i][j]=(float)2.193*buf->Median();
		}
	}
	free(buf);
	cov=covariation;
	return(cov);
}


MatrixObject *ClusterPlus::MeanCovariation(){
int i,j,k;
float x,y;
float *mean;
//MatrixObject *cov;
	cov= new MatrixObject( first->number );
   if(NumberOfElements<2) return(cov);
	mean=(float *)calloc(first->number,sizeof(float));
	if(!mean) return(0);
	for(i=0;i<first->number;i++){
		current=first;
		for(j=0;j<NumberOfElements;j++){
			mean[i]+=current->vec[i];
			current=current->next;
		}
	}
	for(i=0;i<first->number;i++){
		mean[i]=mean[i]/NumberOfElements;
	}
	for(i=0;i<cov->numColumns;i++){
		for(j=0;j<cov->numColumns;j++){
			current=first;
			for(k=0;k<NumberOfElements;k++){
				x=(current->vec[i]-mean[i]);
				y=(current->vec[j]-mean[j]);
				x=x*y;
				cov->matrix[i][j]+=x;
				current=current->next;
			}
			cov->matrix[i][j]/=NumberOfElements-1;
		}
	}
	free(mean);
	return(cov);
}

MatrixObject *ClusterPlus::EstimateCovariation(){
int i,j,k;
float *mean;
MatrixObject *m;
	m= new MatrixObject(first->number);
	mean=(float *)calloc(first->number,sizeof(float));
	if(!mean) return(0);
	for(i=0;i<first->number;i++){
		current=first;
		for(j=0;j<NumberOfElements;j++){
			mean[i]+=current->vec[i];
			current=current->next;
		}
	}
	for(i=0;i<first->number;i++){
		mean[i]=mean[i]/NumberOfElements;
	}
	for(i=0;i<m->numColumns;i++){
		for(j=0;j<m->numColumns;j++){
			current=first;
			for(k=0;k<NumberOfElements;k++){
				m->matrix[i][j]+=current->vec[i]*current->vec[j];
				current=current->next;
			}
			m->matrix[i][j]=m->matrix[i][j]-(first->number*mean[i]*mean[j]);
			m->matrix[i][j]=m->matrix[i][j]/(first->number-1);
		}
	}
	free(mean);
	return(m);
}

MatrixObject *ClusterPlus::Correlation(){
int i,j,k;
float *mean;
float *disp;
MatrixObject *cor;
	cor= new MatrixObject(first->number);
	mean=(float *)calloc(first->number,sizeof(float));
	if(!mean) return(NULL);
	disp=(float *)calloc(first->number,sizeof(float));
	if(!disp) return(NULL);
	for(i=0;i<first->number;i++){
		current=first;
		for(j=0;j<NumberOfElements;j++){
			mean[i]+=current->vec[i];
			current=current->next;
		}
	}
	for(i=0;i<first->number;i++){
		mean[i]=mean[i]/NumberOfElements;
	}

	for(i=0;i<first->number;i++){
		disp[i]=0;
		current=first;
		for(j=0;j<NumberOfElements;j++){
			disp[i]+=(current->vec[i]-mean[i])*(current->vec[i]-mean[i]);
			current=current->next;
		}
	}
	for(i=0;i<first->number;i++){
		disp[i]=(float)sqrt((double)(disp[i]/(NumberOfElements)));
	}
	for(i=0;i<cor->numColumns;i++){
		for(j=0;j<cor->numColumns;j++){
			if(i==j){
				cor->matrix[i][j]=1;
				continue;
			}
			current=first;
			for(k=0;k<NumberOfElements;k++){
				cor->matrix[i][j]+=(current->vec[i]-mean[i])*(current->vec[j]-mean[j]);
				current=current->next;
			}
			cor->matrix[i][j]=cor->matrix[i][j]/NumberOfElements;
			cor->matrix[i][j]=cor->matrix[i][j]/(disp[i]*disp[j]);
		}
	}
	free(mean);
	free(disp);
	return(cor);
}


ClusterPlus *ClusterPlus::PrincipalComponents(int nfirst){
//nfirst - number of first principal components;
int i,j,k,l;
int *remark,*index;
float max;
ClusterPlus *comps;
ClusterElement *ce;
	if(!cov) cov=Covariation();
	if(nfirst==0) nfirst=cov->numRaws;
	remark=(int *)calloc(cov->numRaws,sizeof(int));
	index=(int *)calloc(cov->numRaws,sizeof(int));
   cov->make_eigenvectors();
	comps= new ClusterPlus();

	l=0;
	for(i=0;i<nfirst;i++){
		max=0.;
		for(j=0;j<cov->numRaws;j++){
			if((cov->Eigval->vec[j]>max)&&(!remark[j])){
				max=cov->Eigval->vec[j];
				k=j;
			}
			remark[k]=i+1;
		}
      index[l++]=k;
	}
   free(remark);
	for(i=0;i<NumberOfElements;i++){
		ce = new ClusterElement(nfirst,FLOAT);
		for(j=0;j<nfirst;j++){
      	current=first;
			for(k=0;k<cov->numRaws;k++){
				ce->vec[j]+=cov->Eigvec[index[i]]->vec[k]*current->vec[j];
         }
         current=current->next;
		}
      comps->Add(ce);
   }
	free(index);
	return(comps);
}

float ClusterPlus::MahalanobisTo( ClusterPlus *c, MatrixObject *icov ){
int i,j, flag;
float *d,buf;
float r;
	d=(float *)calloc(first->number,sizeof(float));
   if(!centroid) SetCentroid();
   if(!c->centroid) c->SetCentroid();
	flag=0;
   if(!icov) {
   	if(!cov) Covariation();
      icov=cov->Inverse();
      flag++;
   }
   for(i=0;i<first->number;i++){
   	d[i]=centroid->vec[i]-c->centroid->vec[i];
   }
   r=0; buf=0;
   for(i=0;i<first->number;i++){
		for(j=0;j<first->number;j++){
	   	buf+=d[i]*icov->matrix[i][j];
      }
		for(j=0;j<first->number;j++){
	   	r+=buf*d[j];
      }
   }
   free(d);
	if(flag) {
   	delete icov;
   }
   r=sqrt(fabs(r));
	return(r);
}

float ClusterPlus::KullbackTo( ClusterPlus *c ){
int i,j;
float *d,buf,buf1,buf2;
MatrixObject *icov, *icovc;
	d=(float *)calloc(first->number,sizeof(float));
   if(!centroid) SetCentroid();
   if(!c->centroid) c->SetCentroid();
   if(!cov) Covariation();
   if(!c->cov) c->Covariation();
   for(i=0;i<first->number;i++){
   	d[i]=centroid->vec[i] - c->centroid->vec[i];
   }
   buf=0;
	for(i=0;i<first->number;i++){
   	for(j=0;j<first->number;j++){
      	buf+=d[i]*d[j];//scalar multipication of substructed means
      }
   }
   cov->precision=c->cov->precision=0.0001;
   icov=cov->Inverse();

	icovc=c->cov->Inverse();
   buf1=0;                                            
   for(i=0;i<icov->numRaws;i++){
   	//mult. summ of inverced cov. matrices by scalar mult. of substr. mean vectors
   	buf1+=(icov->matrix[i][i]+icovc->matrix[i][i])*buf;
   }

   //diagonal elements of multiplied substractions of cov and inverced cov matrices
	buf2=0;
   for(i=0;i<cov->numRaws;i++){
     	for(j=0;j<cov->numRaws;j++){
         buf2+=(cov->matrix[i][j]-c->cov->matrix[i][j])*(icov->matrix[j][i]-icovc->matrix[j][i]);
		}
   }
   buf=(buf1+buf2)/2;
   return(buf);
}

void ClusterPlus::WriteText(FILE *f){
int i,j;
	fprintf(f,"\nCluster #%d, seed %s, includes %d element(s),\n",number,name,NumberOfElements);
	fprintf(f,"\nQuality: %f, radius %f, made in %d step(s)\n",quality,radius,steps);
	current=first;
   for(i=0;i<NumberOfElements;i++){
      fprintf(f,"%s\t",current->name);
		for(j=0;j<current->number;j++){
			//fprintf(f,"%f\t",current->vec[j]);
		}
		current=current->next;
		fprintf(f,"\n");
	}
}

void ClusterPlus::WriteFlatText(FILE *f){
int i,j;
	current=first;
   for(i=0;i<NumberOfElements;i++){
      fprintf(f,"%d\t%s\t", number, current->name);
		for(j=0;j<current->number;j++){
			fprintf(f,"%f\t",current->vec[j]);
		}
		current=current->next;
		fprintf(f,"\n");
	}
}

void ClusterPlus::WriteXML(FILE *f){
int i,j;
	fprintf(f,"\n<CLUSTER NUMBER=%d NAME=%s NUM_ELEMENTS=%d QUALITY=%f RADIUS=%f STEPS=%d> \n",number,name,NumberOfElements,quality,radius,steps);
	current=first;
	for(i=0;i<NumberOfElements;i++){
		fprintf(f,"<VECTOR NAME=%s, LENGTH=%d>\n",current->name,first->number);
		for(j=0;j<current->number;j++){
			fprintf(f,"%f\t",current->vec[j]);
		}
		current=current->next;
		fprintf(f,"\n</VECTOR>\n");
	}
	fprintf(f,"</CLUSTER>\n");
}

void ClusterPlus::ReadXML(FILE *f){
int i,n;
float x;
char *sp;
static char ss[1000000], buf[100];
ClusterElement *ce;

	while(1){
	   if(!fgets(ss,1000000,f)) return;
	   if( strstr(ss, "/CLUSTER") ){
	      break;
		}
	   if( strstr(ss, "CLUSTER") ){
	   	sp=strstr(ss,"NUMBER=")+7;
	      number=atoi(sp);
	   	sp=strstr(ss,"NAME=")+5;
	      i=0;
	      while(!isspace(*sp)){
	      	buf[i++]=*sp;
	         sp++;
	      }
			name=(char *)calloc(i+1,sizeof(char));
	      strcpy(name,buf);
	   	sp=strstr(ss,"NUM_ELEMENTS=")+13;
	      n=atoi(sp);
	   	sp=strstr(ss,"QUALITY=")+8;
	      quality=atof(sp);
	   	sp=strstr(ss,"RADIUS=")+7;
	      radius=atof(sp);
	   	sp=strstr(ss,"STEPS=")+6;
	      steps=atoi(sp);
		}
	   if( strstr(ss, "/VECTOR") ) continue;
	   if( strstr(ss, "VECTOR") ){
	   	sp=strstr(ss,"LENGTH=")+7;
	      i=atoi(sp);
			ce= new ClusterElement(i,FLOAT);
	   	sp=strstr(ss,"NAME=")+5;
	      i=0;
	      while(isspace(*sp)){
	      	buf[i++]=*sp;
	         sp++;
	      }
			ce->name=(char *)calloc(i+1,sizeof(char));
	      strcpy(ce->name,buf);
			for(i=0;i<ce->number;i++) {
	   		fscanf(f,"%f",&x);
	   	   ce->vec[i]=x;
		   }
			Add(ce);
      }
   }
   if(NumberOfElements!=n+1) printf("incomplete!");
}




///////////////////////////////////////////////////////////////////////////////////
MatrixObject::MatrixObject(){
	name=NULL;
	number=0;
	numRaws=0;
	numColumns=0;
	matrix=NULL;
	Eigval=NULL;
	Eigvec=NULL;
}

MatrixObject::MatrixObject(int raws){
int i;
	name=NULL;
	number=0;
	Eigval=NULL;
	Eigvec=NULL;
	numRaws=numColumns=raws;
	matrix=(float **)calloc(numRaws, sizeof(float *));
	for(i=0;i<numRaws;i++) matrix[i]=(float *)calloc(numRaws,sizeof(float));
}

MatrixObject::MatrixObject(int raws, int columns){
int i;
	name=NULL;
	number=0;
	numRaws=0;
	numColumns=0;
	Eigval=NULL;
	Eigvec=NULL;
	numRaws=raws;
	numColumns=columns;
	matrix=(float **)calloc(numRaws, sizeof(float *));
	for(i=0;i<numRaws;i++) matrix[i]=(float *)calloc(numColumns,sizeof(float));
}

MatrixObject *MatrixObject::Transpone(){
int i,j;
float x;
	if(!matrix) return NULL;
	for(i=0;i<numRaws;i++){
   	for(j=i+1;j<numColumns;j++){
			x=matrix[i][j];
      	matrix[i][j]=matrix[j][i];
         matrix[j][i]=x;
      }
   }
	return this;
}

MatrixObject *MatrixObject::Inverse(){
int    i,j,m,k,kk;
double q;
double *A;

	A=(double *)calloc(((numRaws*numRaws)/2+numRaws),sizeof(double));
	k=0;
	for(i=0;i<numRaws;i++){
		for(j=0;j<=i;j++){
			A[k]=(double)matrix[j][i];
			k++;
		}
	}

	if (numRaws<1) return(NULL);
	if (numRaws==1){
		if (A[0]==0.0) return(NULL);
	 	A[0]=1./A[0];
      return(this);
	}

     /* factorization of given matrix by upper
	triangular matrix and its transpose   */
	if(A[0]<0.0) return(NULL);
   A[0]=sqrt(A[0]);
	if (A[0]==0.0) return(NULL);
   for (j=2; j<numRaws+1; j++) {
   	k=(j-1)*j/2;
      A[k]=A[k]/A[0];
   }
   for (i=2; i<numRaws+1; i++){
		q=0.0;
      k=(i-1)*i/2;
		for (m=1; m<i; m++) {
      	q+=A[m+k-1]*A[m+k-1];
      }
		kk=k+m;
      q=A[kk-1]-q;
		if (q<0.0) {
      	return(NULL);
      }else{
      	A[kk-1]=sqrt(q);
      }
		if (numRaws!=2){
			if (A[kk-1]==0.0) return(NULL);
			for (j=i+1; j<numRaws+1; j++){
				q=0.0; kk=(j-1)*j/2;
				for (m=1; m<i; m++) {
            	q+=A[m+k-1]*A[m+kk-1];
            }
      		A[i+kk-1]=(A[i+kk-1]-q)/A[i+k-1];
			}
		}
	}

     /* invertion of upper triangular matrix */
   for (i=1; i<numRaws+1; i++){
		k=(i+1)*i/2;
   	if (A[k-1]==0.0) {
      	return(NULL);
      }else{
      	A[k-1]=1./A[k-1];
      }
	}
   for (j=numRaws; j>1; j--){
		k=(j-1)*j/2;
   	for (i=j-1; i>0; i--){
      	q=0.0;
      	for (m=i+1; m<j+1; m++){
      		q+=A[i+(m-1)*m/2-1]*A[m+k-1];
         }
      	kk=(i+1)*i/2;
         A[i+k-1]=-A[kk-1]*q;
      }
	}

   /* calculation  of inverse matrix  by
	multiplicating triangular matrices */
   for (i=1; i<numRaws+1; i++){
   	for (j=i; j<numRaws+1; j++){
   		q=0.0;
	   	for (m=j; m<numRaws+1; m++){
         	k=(m-1)*m/2;
            q+=A[i+k-1]*A[j+k-1];
         }
	 		A[i+(j-1)*j/2-1]=q;
		}
   }
	k=0;
	for(i=0;i<numRaws;i++){
		for(j=0;j<=i;j++){
			matrix[j][i]=matrix[i][j]=(float)A[k];
			k++;
		}
	}
   free(A);
	return this;
}

VectorObject *MatrixObject::Eigenvalues(){
	if(Eigval) return(Eigval);
   make_eigenvectors();
   return(Eigval);
}

VectorObject **MatrixObject::Eigenvectors(){
	if(Eigvec) return(Eigvec);
   make_eigenvectors();
   return(Eigvec);
}

void MatrixObject::make_eigenvectors(){

int    n1,n2,i,j,k,m,is,ii,jj,i1,j1,ij1;
double q0,q,qq,qqq,a0,aii,ajj,aij,aim,ajm,cosfi,sinfi,sinfi2,cosfi2,sin2fi;

double *C,*R;

	C=(double *)calloc((numRaws*numRaws)/2,sizeof(double));
   R=(double *)calloc((numRaws*numRaws)/2,sizeof(double));
//	precision=0.0001;
	k=0;
	for(i=0;i<numRaws;i++){
		for(j=0;j<=i;j++){
			C[k]=(double)matrix[j][i];
			k++;
		}
	}


	if(numRaws<1) return;
	if(numRaws==1){
		C[0]=1.0;
	   free(C);
   	free(R);
      return;
	}

	is=0;
	for(i=1; i<numRaws+1; i++){
		for (j=1; j<numRaws+1; j++){
			is+=1;
			if(i==j){
				C[is-1]=1.0;
			}else{
				C[is-1]=0.0;
			}
		}
	}
	/* beginning of iterations */
	while(1){
	 /* recognition of principal element  (q,ii,jj) */
		ii=1; jj=2; if(R[1]<0.0) q=-R[1]; else q=R[1];
		is=0;
		for (j=1; j<numRaws+1; j++){
			for (i=1; i<j+1; i++){
				is+=1;
				if(i!=j){
					if(R[is-1]<0.0){
						a0=-R[is-1];
					}else{
						a0=R[is-1];
					}
				}
				if(a0>q){
					q=R[is-1]; ii=i; jj=j;
				}
			}
		}

	 /* exit */
		if(q<precision){
			for (i=1; i<numRaws+1; i++){
				R[i-1]=R[i*(i+1)/2-1];
			}
			break;
		}

	 /* calculation of rotation angle  (cosfi, sinfi) */
		i1=ii*(ii+1)/2-1; j1=jj*(jj+1)/2-1;
		aii=R[i1]; ajj=R[j1]; aij=q; q=aii-ajj;
		if(q==0.0){
			cosfi=sqrt(2.)/2.; sinfi=cosfi;
		}else{
			qq=2.*aij; qqq=q*q+qq*qq;
			if(qqq<0.0) {
			   free(C);
			   free(R);
				return;
			}else{
				qqq=sqrt(qqq);
			}
			if(qqq==0.0) {
			   free(C);
			   free(R);
         	return;
         }
			if(q<0.0){
				a0=-q;
			}else{
				a0=q;
			}
			cosfi=(1.+a0/qqq)/2.0;
			if(cosfi<0.0){
			   free(C);
			   free(R);
				return;
			}else{
				cosfi=sqrt(cosfi);
			}
			if(q*qq<0.0){
				q0=-1.;
			}else{
				q0=1.;
			}
			if(qq<0.0){
				a0=-qq;
			}else{
				a0=qq;
			}
			sinfi=2.*cosfi*qqq;
			if (sinfi==0.0) {
			   free(C);
			   free(R);
				return;
			}else{
				sinfi=q0*a0/sinfi;
			}
		}

	 /* rotation for matrix R :  R:=inv(C(i,j))*R*C(i,j) */
		sinfi2=sinfi*sinfi; cosfi2=cosfi*cosfi; sin2fi=2.*sinfi*cosfi;
		R[i1]=aii*cosfi2+aij*sin2fi+ajj*sinfi2;
		ij1=ii+(jj-1)*jj/2-1; R[ij1]=0.0;
		R[j1]=aii*sinfi2-aij*sin2fi+ajj*cosfi2;
		if (numRaws>2){
			for (m=1; m<numRaws+1; m++){
				if(m!=ii) if(m!=jj) {
					if (ii<m) n1=ii; else n1=m;
					if (ii<m) n2=m;  else n2=ii;
					i1=n1+(n2-1)*n2/2-1;
					if (jj<m) n1=jj; else n1=m;
					if (jj<m) n2=m;  else n2=jj;
					j1=n1+(n2-1)*n2/2-1;
					aim=R[i1]; ajm=R[j1];
					R[i1]=aim*cosfi+ajm*sinfi;
					R[j1]=-aim*sinfi+ajm*cosfi;
				}
			}
		}

	 /* rotation for matrix C :  C:=C*C(i,j) */
		ij1=(ii-1)*numRaws-1; is=(jj-1)*numRaws-1;
		for (m=1; m<numRaws+1; m++) {
			i1=m+ij1; j1=m+is;
			aim=C[i1]; ajm=C[j1];
			C[i1]=aim*cosfi+ajm*sinfi;
			C[j1]=-aim*sinfi+ajm*cosfi;
		}
	}/* the end of iterations */

	k=0;
	for(i=0;i<numRaws;i++){
		Eigval->vec[i]=(float)C[i];
		for(j=0;j<numRaws;j++){
			Eigvec[i]->vec[j]=(float)R[k];
			k++;
		}
	}
   free(C);
   free(R);
   return;
}

////////////////////////////////////////////////////////////////////////////////////

ClusterList::ClusterList(ClusterPlus *clu){
	cluster=clu;
	root=last=this;
   previous=next=NULL;
}

ClusterList::~ClusterList(){
	delete(cluster);
}

ClusterList *ClusterList::go_next(){
	if(next){
		return(next);
	}else{
		return(this);
	}
}

ClusterList *ClusterList::go_previous(){
	if(previous) {
		return(previous);
	}else{
		return(this);
	}
}

ClusterList *ClusterList::go_end(){
	return(last);
}

ClusterList *ClusterList::go_root(){
	return(root);
}

ClusterList *ClusterList::append(ClusterPlus *clu){
ClusterList *cl;
ClusterList *c;
	c= cl = new ClusterList(clu);
	cl->previous=last;
	last->next=cl;
	cl->root=root;
	while(1){
		cl->last=c;
		if(!cl->previous) break;
		cl=cl->go_previous();
	}
   return(c);
}

ClusterList *ClusterList::insert(ClusterPlus *clu){
	ClusterList *cl = new ClusterList(clu);
	cl->next=next;
	cl->previous=this;
	next->previous=cl;
	next=cl;
	cl->root=root;
	if(this==last){
		last=cl->last=cl;
	}else{
		cl->last=last;
	}
   return(cl);
}


ClusterList *ClusterList::remove(){
ClusterList *cl;
	if(previous) {
		previous->next=next;
	}else{
		cl=this;
		while(1){
			cl->root=next;
         if(!cl->next) break;
			cl=cl->go_next();
		}
	}
	if(next) {
		next->previous=previous;
	}else{
		cl=root;
		while(1){
			cl->last=previous;;
         if(!cl->next) break;
			cl=cl->go_next();
		}
	}
	return(this);
}

