//Copyright(c) Andrey Ptitsyn, PBRC 2002
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
//#include <values.h>

#define NULL 0
#define MAXFLOAT 999999999999.99

//variable types for data vectors
#define FLOAT 1
#define INT 2
#define BINARY 2

//metric types for inter-object distance
#define HEMMING 1	//Hemming binary distance
#define EUCLID 2	//Euclidean distance
#define CORR 3		//Pearson correlation coefficient
#define CHI2 4		//Chi-square distance
#define ZHIV 5		//Zhivotovsky criterion
#define CULB 6		//Kullback criterion

//estimation type
#define MEAN 0 //0 makes it deafault
#define MEDIAN 1

#define MEAN_DISTANCE_INSIDE 1
#define MEAN_DISTANCE_TOCENTER 2
#define VARIANCE 3
#define STDEV 4
#define STEPS 5


class VectorObject{
public:
	char *name;
	int number;
	float *vec;
	int *ivec;
	int var_type;
	float mean;
	float median;

	VectorObject();
	VectorObject(int, int);
	~VectorObject();

	void Read( FILE * );
	void Write( FILE * );
	float Median();
	float Mean();
   float Min();
   float Max();
   float ScalarMultiply(VectorObject *);
};//end class VectorObject


class MatrixObject{
public:
char *name;
int number;
int numRaws;
int numColumns;
float **matrix;
float precision;
VectorObject *Eigval;
VectorObject **Eigvec;

MatrixObject();
MatrixObject(int rows, int columns);
MatrixObject(int raws);

MatrixObject *Transpone();
MatrixObject *Inverse();
VectorObject *Eigenvalues();
VectorObject **Eigenvectors();


void make_eigenvectors();
};//end class MatrixObject


class ClusterElement:public VectorObject{
public:
	int metric;
	int mark;
	ClusterElement *previous;
	ClusterElement *next;

	ClusterElement(int, int);
	ClusterElement(ClusterElement *ce);
   float DistanceTo( ClusterElement *another );
	int EqualsTo( ClusterElement *);
	void SetDistanceMetric(int metric);

protected:
	int BinaryVariance();
	int BinaryDistanceTo( ClusterElement *);
	float EuclideanDistanceTo( ClusterElement *);
   float CorrelationTo( ClusterElement *);

};

class Cluster{
public:
	int metric;
	char *name;
	int number;
	int NumberOfElements;
	int redundant;
	ClusterElement *current;
	ClusterElement *first;
	ClusterElement *last;

	Cluster();
	void Add( ClusterElement *);
	void Add( Cluster *);
	void Remove( ClusterElement *);
	int Redundant();
	Cluster *Divide();
	void Glue(Cluster *);
	void GetElement( ClusterElement *);
	void PickElementsFrom( Cluster * );
	void MarkAll();
	void UnmarkAll();
	void SetDistanceMetric();
};

class ClusterPlus: public Cluster{
public:
	int estimation;
   int steps;
   int mark;
   float quality;
   float var;
   float radius;

	float MaxDistanceInside();
	float MinDistanceInside();
	ClusterElement *centroid;
	MatrixObject *cov;
	MatrixObject *cor;

	ClusterPlus();
//	ClusterPlus(Cluster *Cluster);
	~ClusterPlus();

	float Variance();
	float stdev();
	float MahalanobisTo( ClusterPlus *, MatrixObject *);
   float KullbackTo(ClusterPlus *);
	float DistanceToPoint(ClusterElement *);
	float DistanceBetween(ClusterPlus *);
	float DistanceInside();
	float Radius();
	MatrixObject *Covariation();
	MatrixObject *Correlation();
	MatrixObject *EstimateCovariation();
	MatrixObject *Eigenvalues();
	VectorObject *Eigenvector();
	ClusterElement *SetCentroid();
	ClusterPlus *PrincipalComponents(int);
	void WriteText(FILE *);
   void WriteFlatText(FILE *);
	void WriteXML(FILE *);
   void ReadXML(FILE *);


	float MeanDistanceInside();
	float MeanDistanceToPoint( ClusterElement *);
	float MeanDistanceBetween( ClusterPlus *);
	float MeanDistanceFromCentroid( ClusterPlus *);
	float MedianDistanceInside();
	float MedianDistanceToPoint( ClusterElement *);
	float MedianDistanceBetween( ClusterPlus *);
	float MedianDistanceFromCentroid(ClusterPlus *);
	ClusterElement *MeanCentroid();
	ClusterElement *MedianCentroid();
	MatrixObject *MedianCovariation();
	MatrixObject *MeanCovariation();
};


class ClusterList{
public:
	int mark;
	ClusterList(ClusterPlus *);
	ClusterList *previous;
	ClusterList *next;
	ClusterList *root;
	ClusterList *last;
	ClusterPlus *cluster;
	~ClusterList();

	ClusterList *go_root();
	ClusterList *go_end();
	ClusterList *go_next();
	ClusterList *go_previous();
	ClusterList *swap(ClusterList *);
	ClusterList *append(ClusterPlus *);
	ClusterList *insert(ClusterPlus *);
	ClusterList *remove();
};


