#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>

#include <vtkDataSetWriter.h>
#include <vtkDoubleArray.h>
#include <cmath>

#define LEFT_IDX 1
#define TIP_IDX 0 
#define RIGHT_IDX 2

#define ZERO_ONLY false 

#define DEBUG_MESSAGES false 
#define DEPOSIT_TRACE true
#define END_TRACE true 
#define COLOR_TRACE true 
#define NEWCORNER_TRACE true
#define SHADING_TRACE true

using std::cerr;
using std::endl;

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         viewDir[0] = 0;
         viewDir[1] = 0;
         viewDir[2] = -1;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 5.3;
         alpha = 7.5;
    };
  

    double lightDir[3]; // The direction of the light source
    double viewDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for diffuse lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;

double ceil441(double f)
{
    return ceil(f-0.00001);
}

double floor441(double f)
{
    return floor(f+0.00001);
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void
WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

typedef enum TriangleType {
	FLAT_BOT = 0,
	FLAT_TOP = 1,
	ARBITRARY = 2,
	LINE = 3
} t_type;

typedef struct coordinates {
	double X;
	double Y;
	double Z;
	double color[3];
	double normal[3];
	double norm;

} vertix;

void vSwap(vertix * a, vertix * b){
	vertix temp;
	size_t size = sizeof temp;
	
	memcpy(&temp, a, size);
	memcpy(a, b, size);
	memcpy(b, &temp, size);

	/*
	temp.X = a->X;
	temp.Y = a->Y;
	temp.Z = a->Z;
	temp.color[0] = a->color[0];
	temp.color[1] = a->color[1];
	temp.color[2] = a->color[2];
	memcpy(a, b, sizeof b);	
	a->X = b->X;
	a->Y = b->Y;
	a->Z = b->Z;
	a->color[0] = b->color[0];
	a->color[1] = b->color[1];
	a->color[2] = b->color[2];
	b->X = temp.X;
	b->Y = temp.Y;
	b->Z = temp.Z;
	b->color[0] = temp.color[0];
	b->color[1] = temp.color[1];
	b->color[2] = temp.color[2];
	*/
}

class Triangle
{
  public:
      double	X[3];
      double	Y[3];
      double	Z[3];
      double	colors[3][3];	
	double	normals[3][3];
      void  	transform();
      vertix*	getVertices();
      double	getLeftSlopeReciprocal();
      double	getRightSlopeReciprocal();
      t_type	getType();
      Triangle();
      Triangle(vertix points[3], t_type newType);
      t_type	type;
      vertix	points_[3];
      bool  	isTransformed;
		
	
      void getColor(double x, double y, double * color, double leftEnd, double rightEnd);
      double getShadingValue(double x, double y, double z, double leftEnd, double rightEnd);
      double getZ(double x, double y, double leftEnd, double rightEnd);
      double getRightLineSlope();
      double getLeftLineSlope();

      double getRightLineOffset(double rightSlope);
      double getLeftLineOffset(double leftSlope);
      
      double getLineEnd(double y, double offset, double slope);
  private:
};

std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = (pt[0]+10)*50.0;
        tris[idx].Y[0] = (pt[1]+10)*50.0;
        tris[idx].Z[0] = (pt[2]-10)*0.05;
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = (pt[0]+10)*50.0;
        tris[idx].Y[1] = (pt[1]+10)*50.0;
        tris[idx].Z[1] = (pt[2]-10)*0.05;
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = (pt[0]+10)*50.0;
        tris[idx].Y[2] = (pt[1]+10)*50.0;
        tris[idx].Z[2] = (pt[2]-10)*0.05;
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 }, 
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 } 
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}



void fillRow(double * zBuffer, unsigned char * imageBuffer, double Y, double startX, double endX, int width, int height, Triangle * tri){
	if (Y >= 0 && Y < height){
		double idx;
		for(idx = ceil441(startX); idx <= floor441(endX); idx++){
			if (idx >= 0 && idx < width){
				double color[3];
				tri->getColor(idx, Y, color, startX, endX);
				double z = tri->getZ(idx, Y, startX, endX);
				double shadingValue = tri->getShadingValue(idx, Y, z, startX, endX);
				
				if (DEBUG_MESSAGES && COLOR_TRACE)
				std::cerr << std::setprecision(10) << "RGB: [" << color[0] << ", " 
				<< color[1] << ", " << color[2]<< "]" << endl;	
				int bufferIdx = 3*(Y*width+idx);
				if (z >= zBuffer[bufferIdx]){
					zBuffer[bufferIdx] = z;
					imageBuffer[bufferIdx] = ceil441(255*std::min(1.00, shadingValue*color[0]));
					imageBuffer[bufferIdx + 1] = ceil441(255*std::min(1.00, shadingValue*color[1]));
					imageBuffer[bufferIdx + 2] = ceil441(255*std::min(1.00, shadingValue*color[2]));
				}
			}
		}
	}
}

Triangle::Triangle(){
      isTransformed = false;
}

double Triangle::getShadingValue(double x, double y, double z, double leftEnd, double rightEnd){
	double tLeft = (((points_[TIP_IDX].Y - points_[LEFT_IDX].Y) == 0) ? 
				0 : (y - points_[LEFT_IDX].Y)/(points_[TIP_IDX].Y - points_[LEFT_IDX].Y));
	double tRight = (((points_[TIP_IDX].Y - points_[RIGHT_IDX].Y) == 0) ? 
				0 : (y - points_[RIGHT_IDX].Y)/(points_[TIP_IDX].Y - points_[RIGHT_IDX].Y));
	if (DEBUG_MESSAGES && SHADING_TRACE)
	std::cerr << "<------------- Interpolating Shading ------------->"
	<< endl
	<< "Target: (" << x << ", " << y << ")" 
	<< endl;


	if (DEBUG_MESSAGES && SHADING_TRACE){
		if (tLeft == 0 || tRight == 0){
			std::cerr << "t: ZERO DETECTED!" << endl;
		} else {
			std::cerr << std::setprecision(16) << "t: [" << tLeft << ", " << tRight << "]" << endl;
		}
	}
	
	if (DEBUG_MESSAGES && SHADING_TRACE)
	std::cerr << std::setprecision(16)
	<< "Original Left Normals: [" << points_[LEFT_IDX].normal[0] << ", " 
	<< points_[LEFT_IDX].normal[1] << ", " << points_[LEFT_IDX].normal[2] << "]" 
	<< endl
	<< "Original Right Normals: [" << points_[RIGHT_IDX].normal[0] << ", " 
	<< points_[RIGHT_IDX].normal[1] << ", " << points_[RIGHT_IDX].normal[2] << "]" 
	<< endl
	<< "Original Top Normals: [" << points_[TIP_IDX].normal[0] << ", " 
	<< points_[TIP_IDX].normal[1] << ", " << points_[TIP_IDX].normal[2] << "]" 
	<< endl
	<< "NORMS: [" << points_[TIP_IDX].norm << ", " 
	<< points_[LEFT_IDX].norm << ", " << points_[RIGHT_IDX].norm << "]" 
	<< endl
	<< "Normalized Left Normals: [" 
	<< points_[LEFT_IDX].normal[0] / points_[LEFT_IDX].norm << ", " 
	<< points_[LEFT_IDX].normal[1] / points_[LEFT_IDX].norm << ", " 
	<< points_[LEFT_IDX].normal[2] /points_[LEFT_IDX].norm  << "]" 
	<< endl
	<< "Normalized Right Normals: [" 
	<< points_[RIGHT_IDX].normal[0] / points_[RIGHT_IDX].norm << ", " 
	<< points_[RIGHT_IDX].normal[1] / points_[RIGHT_IDX].norm << ", " 
	<< points_[RIGHT_IDX].normal[2] /points_[RIGHT_IDX].norm  << "]" 
	<< endl
	<< "Normalized Tip Normals: [" 
	<< points_[TIP_IDX].normal[0] / points_[TIP_IDX].norm << ", " 
	<< points_[TIP_IDX].normal[1] / points_[TIP_IDX].norm << ", " 
	<< points_[TIP_IDX].normal[2] /points_[TIP_IDX].norm  << "]" 
	<< endl;
	

	double leftNormalX = points_[LEFT_IDX].normal[0] 
			+ (tLeft)*(points_[TIP_IDX].normal[0] - points_[LEFT_IDX].normal[0]);
	double leftNormalY = points_[LEFT_IDX].normal[1] 
			+ (tLeft)*(points_[TIP_IDX].normal[1] - points_[LEFT_IDX].normal[1]);
	double leftNormalZ = points_[LEFT_IDX].normal[2] 
			+ (tLeft)*(points_[TIP_IDX].normal[2] - points_[LEFT_IDX].normal[2]);

	double leftNorm = std::sqrt(leftNormalX * leftNormalX + leftNormalY * leftNormalY + leftNormalZ * leftNormalZ);

	leftNormalX /= leftNorm;
	leftNormalY /= leftNorm;
	leftNormalZ /= leftNorm;

	double rightNormalX = points_[RIGHT_IDX].normal[0] 
			+ (tRight)*(points_[TIP_IDX].normal[0]-points_[RIGHT_IDX].normal[0]);
	double rightNormalY = points_[RIGHT_IDX].normal[1] 
			+ (tRight)*(points_[TIP_IDX].normal[1]-points_[RIGHT_IDX].normal[1]);
	double rightNormalZ = points_[RIGHT_IDX].normal[2] 
			+ (tRight)*(points_[TIP_IDX].normal[2]-points_[RIGHT_IDX].normal[2]);

	double rightNorm = std::sqrt(rightNormalX * rightNormalX + rightNormalY * rightNormalY 
						+ rightNormalZ * rightNormalZ);
	
	rightNormalX /= rightNorm;
	rightNormalY /= rightNorm;
	rightNormalZ /= rightNorm;

	double tX = ((rightEnd - leftEnd) == 0) ? 0 : ((x - leftEnd)/(rightEnd - leftEnd));

	double normalX = leftNormalX + (tX)*(rightNormalX - leftNormalX);
	double normalY = leftNormalY + (tX)*(rightNormalY - leftNormalY);
	double normalZ = leftNormalZ + (tX)*(rightNormalZ - leftNormalZ);
	
	double norm = std::sqrt(normalX * normalX + normalY * normalY + normalZ * normalZ); 
	
	normalX /= norm;
	normalY /= norm;
	normalZ /= norm;

	if (DEBUG_MESSAGES && SHADING_TRACE)
	std::cerr << std::setprecision(6) 
	<< "Normals: [" << normalX << ", " << normalY << ", " << normalZ << "]" 
	<< endl
	<< "leftEndNormals: [" << leftNormalX << ", " << leftNormalY << ", " << leftNormalZ << "]" 
	<< endl
	<< "rightEndNormals: [" << rightNormalX << ", " << rightNormalY << ", " << rightNormalZ << "]" << endl;
	
	double leftShadingValue = 0; //points_[LEFT_IDX].Z + (tLeft)*(points_[TIP_IDX].Z - points_[LEFT_IDX].Z);

	double rightShadingValue = 0; //points_[RIGHT_IDX].Z + (tRight)*(points_[TIP_IDX].Z - points_[RIGHT_IDX].Z);
	
	double lDotN = (normalX * lp.lightDir[0] + normalY * lp.lightDir[1] + normalZ * lp.lightDir[2]); 
	
	if (DEBUG_MESSAGES && SHADING_TRACE)
	std::cerr << "LdotN: " 
	<< lDotN 
	<< endl; 
	
	double vectorR[3];
	vectorR[0] = 2 * lDotN * normalX - lp.lightDir[0];	
	vectorR[1] = 2 * lDotN * normalY - lp.lightDir[1];	
	vectorR[2] = 2 * lDotN * normalZ - lp.lightDir[2];	

	double rDotV = (vectorR[0] * lp.viewDir[0] + vectorR[1] * lp.viewDir[1] + vectorR[2] * lp.viewDir[2]);
	
	return lp.Ka + std::abs(lp.Kd * lDotN) + std::max(0.00, lp.Ks * std::pow(rDotV, lp.alpha));	

	//return lp.Ks * std::max(0.00, std::pow(rDotV, lp.alpha));
	//return std::abs(lp.Kd * lDotN);
	//return std::max(0.00, lp.Kd * (normalX * lp.lightDir[0] + normalY * lp.lightDir[1] + normalZ * lp.lightDir[2]));
	//std::max(0.00, (normalX * lp.lightDir[0] + normalY * lp.lightDir[1] + normalZ * lp.lightDir[2]));
	//return lp.Ka + (tX)*(leftShadingValue - rightShadingValue);
	//return 1 + (tX)*(leftShadingValue - rightShadingValue);
}


double Triangle::getZ(double x, double y, double leftEnd, double rightEnd){

	double tLeft = (y - points_[LEFT_IDX].Y)/(points_[TIP_IDX].Y - points_[LEFT_IDX].Y);
	double tRight = (y - points_[RIGHT_IDX].Y)/(points_[TIP_IDX].Y - points_[RIGHT_IDX].Y);

	double leftZ = points_[LEFT_IDX].Z + (tLeft)*(points_[TIP_IDX].Z - points_[LEFT_IDX].Z);

	double rightZ = points_[RIGHT_IDX].Z + (tRight)*(points_[TIP_IDX].Z - points_[RIGHT_IDX].Z);
	
	double tX = (((rightEnd - leftEnd) == 0) ? 0 : (x - leftEnd)/(rightEnd - leftEnd));

	return leftZ + (tX)*(rightZ - leftZ);

}

void Triangle::getColor(double x, double y, double * color, double leftEnd, double rightEnd){

	double tLeft = (((points_[TIP_IDX].Y - points_[LEFT_IDX].Y) == 0) ? 
				0 : (y - points_[LEFT_IDX].Y)/(points_[TIP_IDX].Y - points_[LEFT_IDX].Y));
	double tRight = (((points_[TIP_IDX].Y - points_[RIGHT_IDX].Y) == 0) ? 
				0 : (y - points_[RIGHT_IDX].Y)/(points_[TIP_IDX].Y - points_[RIGHT_IDX].Y));

	double leftR = points_[LEFT_IDX].color[0] + (tLeft)*(points_[TIP_IDX].color[0] - points_[LEFT_IDX].color[0]);
	double leftG = points_[LEFT_IDX].color[1] + (tLeft)*(points_[TIP_IDX].color[1] - points_[LEFT_IDX].color[1]);
	double leftB = points_[LEFT_IDX].color[2] + (tLeft)*(points_[TIP_IDX].color[2] - points_[LEFT_IDX].color[2]);

	double rightR = points_[RIGHT_IDX].color[0] + (tRight)*(points_[TIP_IDX].color[0]-points_[RIGHT_IDX].color[0]);
	double rightG = points_[RIGHT_IDX].color[1] + (tRight)*(points_[TIP_IDX].color[1]-points_[RIGHT_IDX].color[1]);
	double rightB = points_[RIGHT_IDX].color[2] + (tRight)*(points_[TIP_IDX].color[2]-points_[RIGHT_IDX].color[2]);
	
	double tX = ((rightEnd - leftEnd) == 0) ? 0 : ((x - leftEnd)/(rightEnd - leftEnd));

	color[0] = leftR + (tX)*(rightR - leftR);
	color[1] = leftG + (tX)*(rightG - leftG);
	color[2] = leftB + (tX)*(rightB - leftB);

}

Triangle::Triangle(vertix points[3], t_type newType){
	for(int i = 0; i < 3; i++){
		points_[i] = points[i];
	}
	type = newType;
	isTransformed = true;
}

double Triangle::getRightLineSlope(){
	return ((points_[TIP_IDX].X - points_[RIGHT_IDX].X) == 0) ? 0 : (points_[TIP_IDX].Y-points_[RIGHT_IDX].Y)/(points_[TIP_IDX].X-points_[RIGHT_IDX].X);
}

double Triangle::getLeftLineSlope(){
	return ((points_[TIP_IDX].X - points_[LEFT_IDX].X) == 0) ? 0 : (points_[TIP_IDX].Y-points_[LEFT_IDX].Y)/(points_[TIP_IDX].X-points_[LEFT_IDX].X);
}

double Triangle::getRightLineOffset(double rightSlope){
	return ((rightSlope ==0) ? points_[RIGHT_IDX].X : (points_[RIGHT_IDX].Y - rightSlope * points_[RIGHT_IDX].X)); 
}

double Triangle::getLeftLineOffset(double leftSlope){
	return ((leftSlope == 0) ? points_[LEFT_IDX].X : (points_[LEFT_IDX].Y - leftSlope * points_[LEFT_IDX].X)); 
}

double Triangle::getLineEnd(double y, double offset, double slope){
	return ((slope == 0) ? (offset) : (y - offset) / slope);	
}
double Triangle::getLeftSlopeReciprocal(){
	if (!isTransformed){
		transform();
	}
	return ((points_[0].Y - points_[1].Y) == 0) ? 0 : (points_[0].X-points_[1].X)/(points_[0].Y-points_[1].Y);
}

double Triangle::getRightSlopeReciprocal(){
	if (!isTransformed){
		transform();
	}
	return ((points_[0].Y - points_[2].Y) == 0) ? 0 : (points_[0].X-points_[2].X)/(points_[0].Y-points_[2].Y);
}

void Triangle::transform(){
		vertix leftPoint, topPoint, rightPoint;
		int leftIdx, topIdx, rightIdx;
		double maxY = 0;
		
		for (int i = 0; i < 3; i++){
			this->points_[i].X = this->X[i];
			this->points_[i].Y = this->Y[i];
			this->points_[i].Z = this->Z[i];

			for (int j = 0; j < 3; j++){
				this->points_[i].color[j] = this->colors[i][j];
				this->points_[i].normal[j] = this->normals[i][j];
			}
		}
		for (int i = 0; i < 4; i++){
			if (points_[i%2].Y < points_[(i%2)+1].Y){
				vSwap(&points_[i%2], &points_[i%2+1]);
			}
		}
		if (points_[1].Y == points_[2].Y){
			if (points_[1].X > points_[2].X){
				vSwap(&points_[1], &points_[2]);
			}
			type = FLAT_BOT;
		} else if (points_[0].Y == points_[1].Y){
			vSwap(&points_[0], &points_[2]);
			if (points_[1].X > points_[2].X){
				vSwap(&points_[1], &points_[2]);
			}
			type = FLAT_TOP;
		} else {
			type = ARBITRARY;
		}
	
		for (int i = 0; i < 3; i++){
			this->points_[i].norm = std::sqrt((this->points_[i].normal[0] * this->points_[i].normal[0]) 
			+ (this->points_[i].normal[1] * this->points_[i].normal[1]) 
			+ (this->points_[i].normal[2] * this->points_[i].normal[2]));
			this->points_[i].normal[0] /= this->points_[i].norm;
			this->points_[i].normal[1] /= this->points_[i].norm;
			this->points_[i].normal[2] /= this->points_[i].norm;
		}
		
	
		isTransformed = true;
}

t_type Triangle::getType(){
	return type;
}

vertix * Triangle::getVertices(){
	if (!isTransformed){
		transform();
	}

	return points_;
}

class Screen
{
  public:
      unsigned char   *buffer;
      double * zBuffer;
      int width, height;
      void depositTriangle(Triangle triangle);

};

void Screen::depositTriangle(Triangle triangle){
	vertix * temp = triangle.getVertices();
	
	if(DEBUG_MESSAGES && DEPOSIT_TRACE)
   	std::cerr << std::setprecision(10) 
	<< "[[============== depositing new triangle ============]]"
	<< endl << "n[0]: [" << (double)temp[0].normal[0] << ", " 
	<< (double)temp[0].normal[1] << ", " << (double)temp[0].normal[2] << "]"
	<< endl << "n[1]: [" << (double)temp[1].normal[0] << ", " 
	<< (double)temp[1].normal[1] << ", " << (double)temp[1].normal[2] << "]"
	<< endl << "n[2]: [" << (double)temp[2].normal[0] << ", " 
	<< (double)temp[2].normal[1] << ", " << (double)temp[2].normal[2] << "]"
	<< endl << "c[0]: [" << (double)temp[0].color[0] << ", " 
	<< (double)temp[0].color[1] << ", " << (double)temp[0].color[2] << "]"
	<< endl << "c[1]: [" 
	<< (double)temp[1].color[0] << ", " << (double)temp[1].color[1] << ", "
	<< (double)temp[1].color[2] << "]" 
	<< endl << "c[2]: ["
	<< (double)temp[2].color[0] << ", "
	<< (double)temp[2].color[1] << ", " << (double)temp[2].color[2] << "]"
	<< endl << "v[0]: "
	<< "(" << temp[0].X << ",  " << temp[0].Y << ", " << temp[0].Z << ")"
	<< endl << "v[1]: "
	<< "(" << temp[1].X << ",  " << temp[1].Y << ", " << temp[1].Z << ")"
	<< endl << "v[2]: "
	<< "(" << temp[2].X << ",  " << temp[2].Y << ", " << temp[2].Z << ")"
	<< endl
	<< "                ~~~~~~~~~~~~~~~~~~~~~~~              "
	<< endl;
	
	t_type type = triangle.getType();
	double rSlope = triangle.getRightSlopeReciprocal();
	double lSlope = triangle.getLeftSlopeReciprocal();
	
	double rightSlope = triangle.getRightLineSlope();
	double leftSlope = triangle.getLeftLineSlope();

	if (DEBUG_MESSAGES && DEPOSIT_TRACE)
		std::cerr << "Left Slope: " << leftSlope << " Right Slope : " << rightSlope << endl;
	
	double rightOffset = triangle.getRightLineOffset(rightSlope);
	double leftOffset = triangle.getLeftLineOffset(leftSlope);
	
	if (DEBUG_MESSAGES && DEPOSIT_TRACE)
		std::cerr << "Left Offset: " << leftOffset << " Right Offset: " << rightOffset << endl;
		
	if (type == FLAT_BOT){
		for(double idxY = ceil441(temp[LEFT_IDX].Y); idxY <= floor441(temp[TIP_IDX].Y); idxY++){
			double deltaY = idxY - ceil441(temp[LEFT_IDX].Y);
			double leftLineEnd = triangle.getLineEnd(idxY, leftOffset, leftSlope);
			double rightLineEnd = triangle.getLineEnd(idxY, rightOffset, rightSlope);
			
			if (DEBUG_MESSAGES && END_TRACE)
				std::cerr << "Left End: " << leftLineEnd << " Right End: " << rightLineEnd << endl;
		
				
			fillRow(zBuffer, buffer, idxY, leftLineEnd, rightLineEnd, 
				width, height, &triangle);
		}
	} else if (type == FLAT_TOP){
		for(double idxY = floor441(temp[LEFT_IDX].Y); idxY >= ceil441(temp[TIP_IDX].Y); idxY--){
			double deltaY = floor(temp[LEFT_IDX].Y) - idxY;
			double leftLineEnd = triangle.getLineEnd(idxY, leftOffset, leftSlope);
			double rightLineEnd = triangle.getLineEnd(idxY, rightOffset, rightSlope);
			
			if (DEBUG_MESSAGES && END_TRACE)
				std::cerr << "Left End: " << leftLineEnd << " Right End: " << rightLineEnd << endl;
	
			
			fillRow(zBuffer, buffer, idxY, leftLineEnd, rightLineEnd, 
				width, height, &triangle);
		}
	} else if (type == ARBITRARY){
		vertix newCorner;
			
		double slope = triangle.getRightLineSlope();
      		double offset = triangle.getRightLineOffset(slope);
      

		newCorner.Y = triangle.points_[LEFT_IDX].Y;
      		newCorner.X = triangle.getLineEnd(newCorner.Y, offset, slope);

		double newCornerT = (newCorner.Y - triangle.points_[RIGHT_IDX].Y)
				/ (triangle.points_[TIP_IDX].Y - triangle.points_[RIGHT_IDX].Y);
		double newCornerR = triangle.points_[RIGHT_IDX].color[0] 
				+ (newCornerT)*(triangle.points_[TIP_IDX].color[0] 
				- triangle.points_[RIGHT_IDX].color[0]);
		double newCornerG = triangle.points_[RIGHT_IDX].color[1] 
				+ (newCornerT)*(triangle.points_[TIP_IDX].color[1] 
				- triangle.points_[RIGHT_IDX].color[1]);
		double newCornerB = triangle.points_[RIGHT_IDX].color[2] 
				+ (newCornerT)*(triangle.points_[TIP_IDX].color[2] 
				- triangle.points_[RIGHT_IDX].color[2]);
		double newCornerZ = triangle.points_[RIGHT_IDX].Z 
				+ (newCornerT)*(triangle.points_[TIP_IDX].Z 
				- triangle.points_[RIGHT_IDX].Z);
		double newCornerNormalX = triangle.points_[RIGHT_IDX].normal[0] 
				+ (newCornerT)*(triangle.points_[TIP_IDX].normal[0] 
				- triangle.points_[RIGHT_IDX].normal[0]);
		double newCornerNormalY = triangle.points_[RIGHT_IDX].normal[1] 
				+ (newCornerT)*(triangle.points_[TIP_IDX].normal[1] 
				- triangle.points_[RIGHT_IDX].normal[1]);
		double newCornerNormalZ = triangle.points_[RIGHT_IDX].normal[2] 
				+ (newCornerT)*(triangle.points_[TIP_IDX].normal[2] 
				- triangle.points_[RIGHT_IDX].normal[2]);
		
		newCorner.color[0] = newCornerR;
		newCorner.color[1] = newCornerG;
		newCorner.color[2] = newCornerB;
		newCorner.normal[0] = newCornerNormalX;
		newCorner.normal[1] = newCornerNormalY;
		newCorner.normal[2] = newCornerNormalZ;
		newCorner.Z = newCornerZ;

		double norm = std::sqrt(newCornerNormalX * newCornerNormalX 
					+ newCornerNormalY * newCornerNormalY 
					+ newCornerNormalZ * newCornerNormalZ);
		
		newCorner.normal[0] /= norm;
		newCorner.normal[1] /= norm;
		newCorner.normal[2] /= norm;


		if (DEBUG_MESSAGES && NEWCORNER_TRACE)
		std::cerr << "newCorner RGB: [" << newCorner.color[0] << ", " << newCorner.color[1] 
		<< ", " << newCorner.color[2] << "]" << endl;
		
		if (DEBUG_MESSAGES && NEWCORNER_TRACE)
		std::cerr << "newCorner normal: [" << newCorner.normal[0] << ", " << newCorner.normal[1] 
		<< ", " << newCorner.normal[2] << "]" << endl;
		
		//newCorner.Z = triangle.getZ(newCorner.X, newCorner.Y, triangle.points_[LEFT_IDX].X, 
		//				triangle.points_[RIGHT_IDX].X);
		if (DEBUG_MESSAGES && NEWCORNER_TRACE)
		std::cerr << "newCorner Z set to: " << newCorner.Z << endl;
			
		vertix newVertices[3];
		if (triangle.points_[LEFT_IDX].X > newCorner.X){
			memcpy(&newVertices[LEFT_IDX], &newCorner, sizeof newVertices[LEFT_IDX]);
			memcpy(&newVertices[RIGHT_IDX], &triangle.points_[LEFT_IDX], sizeof newVertices[RIGHT_IDX]);
		} else {
			memcpy(&newVertices[RIGHT_IDX], &newCorner, sizeof newVertices[RIGHT_IDX]);
			memcpy(&newVertices[LEFT_IDX], &triangle.points_[LEFT_IDX], sizeof newVertices[LEFT_IDX]);
		}
		memcpy(&newVertices[TIP_IDX], &triangle.points_[RIGHT_IDX], sizeof newVertices[TIP_IDX]);
		Triangle newTriangle(newVertices, FLAT_TOP);
		triangle.type = FLAT_BOT;
		if (triangle.points_[LEFT_IDX].X > newCorner.X){
			memcpy(&triangle.points_[RIGHT_IDX], &newCorner, sizeof triangle.points_[RIGHT_IDX]);
			memcpy(&triangle.points_[LEFT_IDX], &triangle.points_[LEFT_IDX], 
				sizeof triangle.points_[LEFT_IDX]);
		} else {
			memcpy(&triangle.points_[RIGHT_IDX], &newCorner, sizeof triangle.points_[RIGHT_IDX]);
			memcpy(&triangle.points_[LEFT_IDX], &triangle.points_[LEFT_IDX], sizeof newVertices[LEFT_IDX]);
		}
			
		if (triangle.points_[LEFT_IDX].X > triangle.points_[RIGHT_IDX].X){
			vSwap(&triangle.points_[LEFT_IDX], &triangle.points_[RIGHT_IDX]);
		}
		
		if (DEBUG_MESSAGES)
		std::cerr << std::setprecision(10) << "SPLIT INTO:" << endl << "	(" 
		<< triangle.points_[TIP_IDX].X << ", " << triangle.points_[TIP_IDX].Y << "), (" 
		<< triangle.points_[LEFT_IDX].X << ", " << triangle.points_[LEFT_IDX].Y  << "), (" 
		<< triangle.points_[RIGHT_IDX].X << ", " << triangle.points_[RIGHT_IDX].Y << ")" 
		<< endl 
		<< "	(" << newTriangle.points_[TIP_IDX].X << ", " << newTriangle.points_[TIP_IDX].Y 
		<< "), (" << newTriangle.points_[LEFT_IDX].X << ", " << newTriangle.points_[LEFT_IDX].Y  
		<< "), (" << newTriangle.points_[RIGHT_IDX].X << ", " << newTriangle.points_[RIGHT_IDX].Y 
		<< ")" 
		<< endl
		<< "	"
		<< "[ " << triangle.points_[TIP_IDX].color[0] << ", " <<  triangle.points_[TIP_IDX].color[1]
		<< ", " << triangle.points_[TIP_IDX].color[2] << "], [" << triangle.points_[LEFT_IDX].color[0] 
		<< ", " <<  triangle.points_[LEFT_IDX].color[1] << ", " << triangle.points_[LEFT_IDX].color[2]
		<< "], [" << triangle.points_[RIGHT_IDX].color[0] << ", " 
		<< triangle.points_[RIGHT_IDX].color[1] << ", " << triangle.points_[RIGHT_IDX].color[2] << "]" 
		<< endl
		<< "	"
		<< "[ " << newTriangle.points_[TIP_IDX].color[0] << ", " <<  triangle.points_[TIP_IDX].color[1]
                << ", " << newTriangle.points_[TIP_IDX].color[2] << "], [" 
		<< triangle.points_[LEFT_IDX].color[0] << ", " <<  newTriangle.points_[LEFT_IDX].color[1] 
		<< ", " << triangle.points_[LEFT_IDX].color[2] << "], [" 
		<< newTriangle.points_[RIGHT_IDX].color[0] << ", " << newTriangle.points_[RIGHT_IDX].color[1] 
		<< ", " << triangle.points_[RIGHT_IDX].color[2] << "]"
		<< endl; 
		
		depositTriangle(triangle);
		depositTriangle(newTriangle);
		
	}

}

int main()
{
   int imgWidth = 1000, imgHeight = 1000;
   vtkImageData *image = NewImage(imgWidth, imgHeight);
   int imgSize= imgWidth*imgHeight;
   int bufferSize = 3*imgSize; 
   
   unsigned char *buffer = 
     (unsigned char *) image->GetScalarPointer(0,0,0);
   for (int i = 0 ; i < bufferSize; i++)
       buffer[i] = 0;

  
   double * zBuffer = (double *) malloc(bufferSize * sizeof(double));
   for (int i = 0 ; i < bufferSize; i++)
       zBuffer[i] = -1;
  

   std::vector<Triangle> triangles = GetTriangles();
   
   Screen screen;
   screen.buffer = buffer;
   screen.width = imgWidth;
   screen.height = imgHeight;
   screen.zBuffer = zBuffer;
	
   std::cerr << "Depositing";
   Triangle test;
   test.X[0] = 600;
   test.X[1] = 650;
   test.X[2] = 550;
   test.Y[0] = 600;
   test.Y[1] = 550;
   test.Y[2] = 500;
   test.Z[0] = 1;
   test.Z[1] = 1;
   test.Z[2] = 1;
   test.colors[0][0] = 1;
   test.colors[0][1] = 0;
   test.colors[0][2] = 0;
   test.colors[1][0] = 0;
   test.colors[1][1] = 1;
   test.colors[1][2] = 0;
   test.colors[2][0] = 0;
   test.colors[2][1] = 0;
   test.colors[2][2] = 1;
   Triangle test2;
   test2.X[0] = 500;
   test2.X[1] = 450;
   test2.X[2] = 460;
   test2.Y[0] = 600;
   test2.Y[1] = 550;
   test2.Y[2] = 500;
   test2.Z[0] = 1;
   test2.Z[1] = 1;
   test2.Z[2] = 1;
   test2.colors[0][0] = 1;
   test2.colors[0][1] = 0;
   test2.colors[0][2] = 0;
   test2.colors[1][0] = 0;
   test2.colors[1][1] = 1;
   test2.colors[1][2] = 0;
   test2.colors[2][0] = 0;
   test2.colors[2][1] = 0;
   test2.colors[2][2] = 1;
   Triangle test3;
   test3.X[0] = 600;
   test3.X[1] = 650;
   test3.X[2] = 550;
   test3.Y[0] = 450;
   test3.Y[1] = 400;
   test3.Y[2] = 400;
   test3.Z[0] = 1;
   test3.Z[1] = 1;
   test3.Z[2] = 1;
   test3.colors[0][0] = 1;
   test3.colors[0][1] = 0;
   test3.colors[0][2] = 0;
   test3.colors[1][0] = 0;
   test3.colors[1][1] = 1;
   test3.colors[1][2] = 0;
   test3.colors[2][0] = 0;
   test3.colors[2][1] = 0;
   test3.colors[2][2] = 1;
   Triangle test4;
   test4.X[0] = 450;
   test4.X[1] = 500;
   test4.X[2] = 400;
   test4.Y[0] = 400;
   test4.Y[1] = 450;
   test4.Y[2] = 450;
   test4.Z[0] = 1;
   test4.Z[1] = 1;
   test4.Z[2] = 1;
   test4.colors[0][0] = 1;
   test4.colors[0][1] = 0;
   test4.colors[0][2] = 0;
   test4.colors[1][0] = 0;
   test4.colors[1][1] = 1;
   test4.colors[1][2] = 0;
   test4.colors[2][0] = 0;
   test4.colors[2][1] = 0;
   test4.colors[2][2] = 1;
   Triangle test5;
   test5.X[0] = 300;
   test5.X[1] = 300;
   test5.X[2] = 250;
   test5.Y[0] = 400;
   test5.Y[1] = 450;
   test5.Y[2] = 450;
   test5.Z[0] = 1;
   test5.Z[1] = 1;
   test5.Z[2] = 1;
   test5.colors[0][0] = 1;
   test5.colors[0][1] = 0;
   test5.colors[0][2] = 0;
   test5.colors[1][0] = 0;
   test5.colors[1][1] = 1;
   test5.colors[1][2] = 0;
   test5.colors[2][0] = 0;
   test5.colors[2][1] = 0;
   test5.colors[2][2] = 1;
   Triangle test6;
   test6.X[0] = 300;
   test6.X[1] = 350;
   test6.X[2] = 300;
   test6.Y[0] = 550;
   test6.Y[1] = 600;
   test6.Y[2] = 600;
   test6.Z[0] = 1;
   test6.Z[1] = 1;
   test6.Z[2] = 1;
   test6.colors[0][0] = 1;
   test6.colors[0][1] = 0;
   test6.colors[0][2] = 0;
   test6.colors[1][0] = 0;
   test6.colors[1][1] = 1;
   test6.colors[1][2] = 0;
   test6.colors[2][0] = 0;
   test6.colors[2][1] = 0;
   test6.colors[2][2] = 1;
   Triangle test7;
   test7.X[0] = 150;
   test7.X[1] = 200;
   test7.X[2] = 150;
   test7.Y[0] = 550;
   test7.Y[1] = 550;
   test7.Y[2] = 600;
   test7.Z[0] = 1;
   test7.Z[1] = 1;
   test7.Z[2] = 1;
   test7.colors[0][0] = 1;
   test7.colors[0][1] = 0;
   test7.colors[0][2] = 0;
   test7.colors[1][0] = 0;
   test7.colors[1][1] = 1;
   test7.colors[1][2] = 0;
   test7.colors[2][0] = 0;
   test7.colors[2][1] = 0;
   test7.colors[2][2] = 1;
   Triangle test8;
   test8.X[0] = 150;
   test8.X[1] = 200;
   test8.X[2] = 200;
   test8.Y[0] = 400;
   test8.Y[1] = 400;
   test8.Y[2] = 450;
   test8.Z[0] = 1;
   test8.Z[1] = 1;
   test8.Z[2] = 1;
   test8.colors[0][0] = 1;
   test8.colors[0][1] = 0;
   test8.colors[0][2] = 0;
   test8.colors[1][0] = 0;
   test8.colors[1][1] = 1;
   test8.colors[1][2] = 0;
   test8.colors[2][0] = 0;
   test8.colors[2][1] = 0;
   test8.colors[2][2] = 1;
   Triangle test9;
   test9.X[0] = 500;
   test9.X[1] = 300;
   test9.X[2] = 400;
   test9.Y[0] = 700;
   test9.Y[1] = 750;
   test9.Y[2] = 750;
   test9.Z[0] = 1;
   test9.Z[1] = 1;
   test9.Z[2] = 1;
   test9.colors[0][0] = 1;
   test9.colors[0][1] = 0;
   test9.colors[0][2] = 0;
   test9.colors[1][0] = 0;
   test9.colors[1][1] = 1;
   test9.colors[1][2] = 0;
   test9.colors[2][0] = 0;
   test9.colors[2][1] = 0;
   test9.colors[2][2] = 1;
   Triangle test10;
   test10.X[0] = 50;
   test10.X[1] = 100;
   test10.X[2] = 200;
   test10.Y[0] = 700;
   test10.Y[1] = 750;
   test10.Y[2] = 750;
   test10.Z[0] = 1;
   test10.Z[1] = 1;
   test10.Z[2] = 1;
   test10.colors[0][0] = 1;
   test10.colors[0][1] = 0;
   test10.colors[0][2] = 0;
   test10.colors[1][0] = 0;
   test10.colors[1][1] = 1;
   test10.colors[1][2] = 0;
   test10.colors[2][0] = 0;
   test10.colors[2][1] = 0;
   test10.colors[2][2] = 1;
   Triangle test11;
   test11.X[0] = 550;
   test11.X[1] = 600;
   test11.X[2] = 650;
   test11.Y[0] = 700;
   test11.Y[1] = 700;
   test11.Y[2] = 750;
   test11.Z[0] = 1;
   test11.Z[1] = 1;
   test11.Z[2] = 1;
   test11.colors[0][0] = 1;
   test11.colors[0][1] = 0;
   test11.colors[0][2] = 0;
   test11.colors[1][0] = 0;
   test11.colors[1][1] = 1;
   test11.colors[1][2] = 0;
   test11.colors[2][0] = 0;
   test11.colors[2][1] = 0;
   test11.colors[2][2] = 1;
   Triangle test12;
   test12.X[0] = 750;
   test12.X[1] = 800;
   test12.X[2] = 700;
   test12.Y[0] = 700;
   test12.Y[1] = 700;
   test12.Y[2] = 750;
   test12.Z[0] = 1;
   test12.Z[1] = 1;
   test12.Z[2] = 1;
   test12.colors[0][0] = 1;
   test12.colors[0][1] = 0;
   test12.colors[0][2] = 0;
   test12.colors[1][0] = 0;
   test12.colors[1][1] = 1;
   test12.colors[1][2] = 0;
   test12.colors[2][0] = 0;
   test12.colors[2][1] = 0;
   test12.colors[2][2] = 1;
   Triangle test13;
   test13.X[0] = 800;
   test13.X[1] = 800;
   test13.X[2] = 700;
   test13.Y[0] = 400;
   test13.Y[1] = 450;
   test13.Y[2] = 350;
   test13.Z[0] = 1;
   test13.Z[1] = 1;
   test13.Z[2] = 1;
   test13.colors[0][0] = 1;
   test13.colors[0][1] = 0;
   test13.colors[0][2] = 0;
   test13.colors[1][0] = 0;
   test13.colors[1][1] = 1;
   test13.colors[1][2] = 0;
   test13.colors[2][0] = 0;
   test13.colors[2][1] = 0;
   test13.colors[2][2] = 1;
   Triangle test14;
   test14.X[0] = 800;
   test14.X[1] = 800;
   test14.X[2] = 700;
   test14.Y[0] = 500;
   test14.Y[1] = 550;
   test14.Y[2] = 600;
   test14.Z[0] = 1;
   test14.Z[1] = 1;
   test14.Z[2] = 1;
   test14.colors[0][0] = 1;
   test14.colors[0][1] = 0;
   test14.colors[0][2] = 0;
   test14.colors[1][0] = 0;
   test14.colors[1][1] = 1;
   test14.colors[1][2] = 0;
   test14.colors[2][0] = 0;
   test14.colors[2][1] = 0;
   test14.colors[2][2] = 1;
   Triangle test15;
   test15.X[0] = 600;
   test15.X[1] = 500;
   test15.X[2] = 500;
   test15.Y[0] = 400;
   test15.Y[1] = 300;
   test15.Y[2] = 250;
   test15.Z[0] = 1;
   test15.Z[1] = 1;
   test15.Z[2] = 1;
   test15.colors[0][0] = 1;
   test15.colors[0][1] = 0;
   test15.colors[0][2] = 0;
   test15.colors[1][0] = 0;
   test15.colors[1][1] = 1;
   test15.colors[1][2] = 0;
   test15.colors[2][0] = 0;
   test15.colors[2][1] = 0;
   test15.colors[2][2] = 1;
   Triangle test16;
   test16.X[0] = 400;
   test16.X[1] = 300;
   test16.X[2] = 300;
   test16.Y[0] = 250;
   test16.Y[1] = 350;
   test16.Y[2] = 300;
   test16.Z[0] = 1;
   test16.Z[1] = 1;
   test16.Z[2] = 1;
   test16.colors[0][0] = 1;
   test16.colors[0][1] = 0;
   test16.colors[0][2] = 0;
   test16.colors[1][0] = 0;
   test16.colors[1][1] = 1;
   test16.colors[1][2] = 0;
   test16.colors[2][0] = 0;
   test16.colors[2][1] = 0;
   test16.colors[2][2] = 1;
/*   screen.depositTriangle(test);
   screen.depositTriangle(test2);
   screen.depositTriangle(test3);
   screen.depositTriangle(test4);
   screen.depositTriangle(test5);
   screen.depositTriangle(test6);
   screen.depositTriangle(test7);
   screen.depositTriangle(test8);
   screen.depositTriangle(test9);
   screen.depositTriangle(test10);
   screen.depositTriangle(test11);
   screen.depositTriangle(test12);
   screen.depositTriangle(test13);
   screen.depositTriangle(test14);
   screen.depositTriangle(test15);
   screen.depositTriangle(test16);
  */ int ID = 4533;
   //std::cerr << std::setprecision(10) << "New: " << (int)triangles[ID].colors[0][0] << ", " << (int)triangles[ID].colors[0][1] << ", " << (int)triangles[ID].colors[0][2] << ", " << (int)triangles[ID].colors[1][0] << ", " << (int)triangles[ID].colors[1][1] << ", " << (int)triangles[ID].colors[1][2] << ", " << (int)triangles[ID].colors[2][0] << ", " << (int)triangles[ID].colors[2][1] << ", " << (int)triangles[ID].colors[2][2] << endl; 
   int clock = (triangles.size()-1)/10000;
   int counter = 1;
   if (ZERO_ONLY){
  	screen.depositTriangle(triangles[0]);
   } else {
   	for (int i = 0; i <= triangles.size() - 1; i++){
   		screen.depositTriangle(triangles[i]);
		if (i%clock == counter){ 
			std::cerr << ".";
			counter++;
		}
	}
   }
   
   std:cerr << endl;
   std::cerr << "Finished Depositing!" << endl;
   WriteImage(image, "output");
}

