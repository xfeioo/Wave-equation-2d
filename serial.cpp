#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <omp.h>
using namespace std;


class Array2D{
private:
    double **tb;
    int row,col;

public:
    Array2D(int r,int c);
    ~Array2D();
    int getRow() const;
    int getCol() const;
    double **getTb() const;
};

class Grids{
private:
    double dx; /**< dx of grids */
    double dy;
    int size[2]; /**< size of grids */
    Array2D *gridNew;
    Array2D *gridOld;
    Array2D *grid;

public:
    Grids(int imax,int jmax,double dx,double dy);
    ~Grids();

    void init(double x,double y,double r,double h);
    void update(double dt,double c);
    void output(int t);
    void applyBondary();

};




int main(int argc,char **argv)
{

    //网格数
    int gridNum[2]={400,400};
    //网格数
    double width=10,length=10;
    double ts=0;
    double te=15;
    //是否输出结果文件
    bool out=false;
    //输出结果文件间隔
    int outFreq = 10;

    double dx=width/gridNum[0];
    double dy=length/gridNum[1];
    double mindx = dx<dy?dx:dy;
    double c=1;
    double dt=0.1*mindx/c;

    double nt=int((te-ts)/dt);
    Grids grid(gridNum[0],gridNum[1],dx,dy);
    //初始化
    grid.init(5,5,1,1);

    double tic,toc;
    //开始计时
    tic=omp_get_wtime();

    for(int t=0;t<nt;t++)
    {
        grid.applyBondary();
        grid.update(dt,c);

        if(t%outFreq==0 and out)
            grid.output(t);

        std::cout<<"Step "<<t<<" at time "<<t*dt<<std::endl;
    }
    //结束计时
    toc=omp_get_wtime();
    std::cout<<"Time use "<<toc-tic<<" s"<<std::endl;

    return 0;
}

int Array2D::getRow() const
{
    return row;
}

int Array2D::getCol() const
{
    return col;
}

double **Array2D::getTb() const
{
    return tb;
}

Array2D::Array2D(int r, int c)
{
    tb=new double*[r];
    double *arr = new double[r*c];
    memset(arr,0,sizeof(double)*r*c);
    for(int i=0;i<r;i++)
    {
        tb[i]=&arr[i*c];
    }
}

Array2D::~Array2D()
{
    delete [] &tb[0][0];
    delete [] &tb[0];
}

Grids::Grids(int imax, int jmax, double dx, double dy)
{
    this->dx=dx;
    this->dy=dy;
    size[0]=imax;
    size[1]=jmax;

    grid = new Array2D(imax,jmax);
    gridNew = new Array2D(imax,jmax);
    gridOld = new Array2D(imax,jmax);
}

Grids::~Grids()
{
    delete grid;
    delete gridNew;
    delete gridOld;

}

void Grids::init(double x, double y, double r, double h)
{
    double **tb = grid->getTb();
    double coord[2];
    const double pi=3.1415926535;
    for(int i=0;i<size[0];i++)
    {
        for(int j=0;j<size[1];j++)
        {
            coord[0]=i*dx;
            coord[1]=j*dy;

            double dist = sqrt((coord[0]-x)*(coord[0]-x)+(coord[1]-y)*(coord[1]-y));
            if(dist<r)
            {
                tb[i][j]=(cos(pi/r*dist)+1)*h;
            }
        }
    }
}

void Grids::update(double dt, double c)
{
    double v;
    double **tb = grid->getTb();
    double **tbOld = gridOld->getTb();
    double **tbNew = gridNew->getTb();

    for(int i=1;i<size[0]-1;i++)
    {
        for(int j=1;j<size[1]-1;j++)
        {
            v = dt*dt*c*c * (
                        (tb[i+1][j]-2.0*tb[i][j]+tb[i-1][j])/dx/dx
                        + (tb[i][j+1]-2.0*tb[i][j]+tb[i][j-1])/dy/dy
                        )
                    +2.0*tb[i][j]-tbOld[i][j];
            tbNew[i][j]=v;
        }
    }

    for(int i=1;i<size[0]-1;i++)
    {
        for(int j=1;j<size[1]-1;j++)
        {
            tbOld[i][j]=tb[i][j];
            tb[i][j]=tbNew[i][j];
        }
    }

}

void Grids::output(int t)
{

    double **tb = grid->getTb();
    string filename="./result/"+std::to_string(t)+".dat";
    FILE *fp=fopen(filename.c_str(),"wb");
    fwrite(&tb[0][0],sizeof(double),size[0]*size[1],fp);
    fclose(fp);
}

void Grids::applyBondary()
{

}
