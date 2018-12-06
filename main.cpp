#include <stdio.h>
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <glut.h>
#include <vector>
#include <set>
#include <algorithm>
using namespace std;

#define pi 3.1415926
#define threHold 0.001
const int MAXN=550;  
const double eps=1e-8; 





/*---------------------我---是---分---割---线------------（结构体定义）------------*/

//空间上任何一个点信息
struct Point {
	double x, y, z;
	Point(){} 
	Point(double xx,double yy,double zz):x(xx),y(yy),z(zz){}
	//两向量之差  
    Point operator -(const Point p1)  
    {  
        return Point(x-p1.x,y-p1.y,z-p1.z);  
    } 

	//两向量之和  
    Point operator +(const Point p1)  
    {  
        return Point(x+p1.x,y+p1.y,z+p1.z);  
    }  

    //叉乘  
    Point operator *(const Point p)  
    {  
        return Point(y*p.z-z*p.y,z*p.x-x*p.z,x*p.y-y*p.x);  
    }  

    Point operator *(double d)  
    {  
        return Point(x*d,y*d,z*d);  
    }  

    Point operator / (double d)  
    {  
        return Point(x/d,y/d,z/d);  
    }  

    //点乘  
    double  operator ^(Point p)  
    {  
        return (x*p.x+y*p.y+z*p.z);  
    } 
};

struct CH3D  
{  
    struct face  
    {  
        //表示凸包一个面上的三个点的编号  
        int a,b,c;  
        //表示该面是否属于最终凸包上的面  
        bool ok;  
    };  

    //初始顶点数  
    int n;  

    //初始顶点  
    Point P[MAXN];  

    //凸包表面的三角形数  
    int num;  

    //凸包表面的三角形  
    face F[8*MAXN];  

    //凸包表面的三角形  
	//g[i][j]存储的是第
    int g[MAXN][MAXN]; 

	//共面点集合，一维是集合数，二维是共面的点数
	vector<set<int>> count;
	vector<Point> polygons[MAXN][2];

    //向量长度  
    double vlen(Point a)  
    {  
        return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);  
    }  

    //叉乘  
    Point cross(const Point &a,const Point &b,const Point &c)  
    {  
        return Point((b.y-a.y)*(c.z-a.z)-(b.z-a.z)*(c.y-a.y),  
                     (b.z-a.z)*(c.x-a.x)-(b.x-a.x)*(c.z-a.z),  
                     (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x)  
                     );  
    }  

    //三角形面积*2  
    double area(Point a,Point b,Point c)  
    {  
        return vlen((b-a)*(c-a));  
    }  

    //四面体有向体积*6  
    double volume(Point a,Point b,Point c,Point d)  
    {  
        return (b-a)*(c-a)^(d-a);  
    }  

    //正：点在面同向  返回的是点和面三点构成的体积，可正也可负
    double dblcmp(Point &p,face &f)  
    {  
        Point m=P[f.b]-P[f.a];  
        Point n=P[f.c]-P[f.a];  
        Point t=p-P[f.a];  
        return (m*n)^t;  
    }  

    void deal(int p,int a,int b)  
    {  
        int f=g[a][b];//搜索与该边相邻的另一个平面  
        face add;  
        if(F[f].ok)  
        {  
            if(dblcmp(P[p],F[f])>eps)  
              dfs(p,f);  
            else  
            {  
                add.a=p;  
                add.b=b;  
                add.c=a;//这里注意顺序，要成右手系  
                add.ok=true;  
                g[p][b]=g[a][p]=g[b][a]=num;  
                F[num++]=add;  
            }  
        }  
    }  

    void dfs(int p,int now)//递归搜索所有应该从凸包内删除的面  
    {  
         F[now].ok=0;  
         deal(p,F[now].b,F[now].a);  
         deal(p,F[now].c,F[now].b);  
         deal(p,F[now].a,F[now].c);  
    }  

	//F[s]和F[t]是否是共面；
    bool same(int s,int t)  
    {  
        Point &a=P[F[s].a];  
        Point &b=P[F[s].b];  
        Point &c=P[F[s].c];  
        return fabs(volume(a,b,c,P[F[t].a]))<eps &&  
               fabs(volume(a,b,c,P[F[t].b]))<eps &&  
               fabs(volume(a,b,c,P[F[t].c]))<eps;  
    }  

    //构建三维凸包  
    void create()  
    {  
        int i,j,tmp;  
        face add;  

        num=0;  
        if(n<4)return;  
    //**********************************************  
    //此段是为了保证前四个点不共面  
        bool flag=true;  
        for(i=1;i<n;i++)  
        {  
            if(vlen(P[0]-P[i])>eps)  
            {  
                swap(P[1],P[i]);  
                flag=false;  
                break;  
            }  
        }  
        if(flag) return;  //所有点都和P[0]点重合
        flag=true;  
        //使前三个点不共线  
        for(i=2;i<n;i++)  
        {  
            if(vlen((P[0]-P[1])*(P[1]-P[i]))>eps)  
            {  
                swap(P[2],P[i]);  
                flag=false;  
                break;  
            }  
        }  
        if(flag)return;  
        flag=true;  
        //使前四个点不共面  
        for(int i=3;i<n;i++)  
        {  
            if(fabs((P[0]-P[1])*(P[1]-P[2])^(P[0]-P[i]))>eps)  
            {  
                swap(P[3],P[i]);  
                flag=false;  
                break;  
            }  
        }  
        if(flag)return;  

    //*****************************************  
        for(i=0;i<4;i++)  
        {  
            add.a=(i+1)%4;  
            add.b=(i+2)%4;  
            add.c=(i+3)%4;  
            add.ok=true;  
            if(dblcmp(P[i],add)>0)
				swap(add.b,add.c);

            g[add.a][add.b]=g[add.b][add.c]=g[add.c][add.a]=num;  
            F[num++]=add;  
        }  
        for(i=4;i<n;i++)  
        {  
            for(j=0;j<num;j++)  
            {  
                if(F[j].ok && dblcmp(P[i],F[j])>eps)  
                {  
                    dfs(i,j);  
                    break;  
                }  
            }  
        }  
        tmp=num;  
        for(i=num=0;i<tmp;i++)  
          if(F[i].ok)  
            F[num++]=F[i];  
		polygon();
    }  

    //表面积  
    double area()  
    {  
        double res=0;  
        if(n==3)  
        {  
            Point p=cross(P[0],P[1],P[2]);  
            res=vlen(p)/2.0;  
            return res;  
        }  
        for(int i=0;i<num;i++)  
          res+=area(P[F[i].a],P[F[i].b],P[F[i].c]);  
        return res/2.0;  
    }  

	//计算凸多面体体积
    double volume()  
    {  
        double res=0;  
        Point tmp(0,0,0);  
        for(int i=0;i<num;i++)  
           res+=volume(tmp,P[F[i].a],P[F[i].b],P[F[i].c]);  
        return fabs(res/6.0);  
    }  

    //表面三角形个数  
    int triangle()  
    {  
        return num;  
    }  

    //表面多边形个数  
    int polygon()  
    {  
        int i,j,res,flag;  
		vector<int> index; //index[i]是第i个三角形是所在的count二维数组下标
		index.resize(num);
        for(i=res=0;i<num;i++)
		{

			flag = 1;
			for(j=0; j<i; j++) {
				if(same(i, j)) {
					count[index[j]].insert(F[i].a);
					count[index[j]].insert(F[i].b);
					count[index[j]].insert(F[i].c);
					index[i] = index[j];
					flag = 0;
					break;
				}
			}
			if(flag) {
				set<int> tempSet;
				tempSet.insert(F[i].a);
				tempSet.insert(F[i].b);
				tempSet.insert(F[i].c);
				count.push_back(tempSet);
				index[i] = count.size()-1;
			}
        }  
		return count.size();
    } 

    //三维凸包重心  
    Point barycenter()  
    {  
        Point ans(0,0,0),o(0,0,0);  
        double all=0;  
        for(int i=0;i<num;i++)  
        {  
            double vol=volume(o,P[F[i].a],P[F[i].b],P[F[i].c]);  
            ans=ans+(o+P[F[i].a]+P[F[i].b]+P[F[i].c])/4.0*vol;  
            all+=vol;  
        }  
        ans=ans/all;  
        return ans;  
    }  

    //点到面的距离  
    double ptoface(Point p,int i)  
    {  
        return fabs(volume(P[F[i].a],P[F[i].b],P[F[i].c],p)/vlen((P[F[i].b]-P[F[i].a])*(P[F[i].c]-P[F[i].a])));  
    }  
}; 

//每个壳体平面
struct surf {
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double a,b,c,d;
	bool isGate;
	double gxmin, gxmax, gymin, gymax, gzmin, gzmax;
	vector<Point> allPoints;
	vector<Point> polygon[2];
};

//每条折射路径
struct lines {
	vector<Point> inner;
	vector<int> surfIndex;
};






/*---------------------我---是---分---割---线------------（函数声明）------------*/


/*glut辅助函数*/
void init(void); //初始化函数
void reshape(int w, int h); //重绘制函数
void keyboard(unsigned char key, int x, int y); //键盘设置函数
void mouse(int button, int state, int x, int y);
void motion(int x, int y); //捕捉鼠标移动函数


/*功能函数*/
void selectPoints(); // 从所有出射光线里筛选从窗口入射的直线函数（如果一开始出射的光线是球面比如）
bool cmp(Point a, Point b);//空间两个点比较函数
void drawPoint(Point p); //绘制一个点函数
void drawLines(Point p1, Point p2); //绘制两点直线函数
double DistanceOfPointToLine(Point* a, Point* b, Point* s); // 计算空间s点到直线ab的距离；
double distanceOfTwoPoints(Point a, Point b); // 计算空间两个点a，b的空间距离；
void drawTriangle(Point a, Point b, Point c); //绘制三角面片的功能函数


/*数据准备部分*/
void dataPrepare(); //数据准备函数
void initSurfs(); //初始化环境平面信息；
void initSingleLight(); //产生单根入射线的方法
void initMultipleLight(); //测试入射线的函数
void analyze(); //循环折射函数
void explicite(); //对在平面内的折射点进行包络体求解的算法
void sortAndConstructThePolygonPoints(); //对包络平面的点进行进一步处理函数
void polygons(); //


/*循环折射计算部分*/
int collision(Point p1, Point p2, Point *hp, Point *next, int index); //碰撞检验函数
Point getCross(Point p1, Point p2, surf s); //计算p1和p2与面域s之间的交点坐标函数，返回交点函数
int validate(Point cross, surf temp); //验证交点是否在某平面之内的函数


/*包络体计算部分*/
void dealWith(vector<Point> &allPoints, vector<Point> polygon[2], Point n1); //求平面内的最大包络多边形
void FindPoint2(vector<Point> &p, Point a, Point b, Point mid, vector<Point> &polygon, Point &n1); //平面求包主题算法



/*下面是逻辑函数*/
void display(void); //主显示函数


/*环境情况绘制*/
void drawEnvironments(); //绘制环境信息
void drawCoordinates(void); //绘制坐标轴函数
void drawEnvironmentLines(); //绘制立方体边框函数
void drawSurfaces(); // 绘制环境平面信息

/*折射情况绘制*/
void drawRefections(); //后期绘制折射所有信息
void drawRefectionTracks(); //绘制所有入射光线的折射路径
void drawRefectionPoints(); //绘制所有光线的所有折射点
void drawRefectionPolygonLines(); //仅绘制在平面上的包络体的最外形轮廓
void drawReflectionPointsInSurfs(); //绘制所有折射路径和所有平面的有效交点
void drawPolygonSurfs(); //绘制在内壳上的最大包络平面（用红色高亮）
void drawSpacePolygonFaces(); //绘制空间凸包络体





/*---------------------我---是---分---割---线--------------（变量定义和初始化）----------*/

GLint mx,my; //position of mouse
GLint m_state=0; //mouse usage
GLfloat x_angle=0.0f, y_angle=0.0f; //angle of eye
GLfloat dist=10.0f; //distance from the eye
double lineWidth = 0.01, pointSize = 0.1;
vector<Point> vp, vtemp[2];
vector<surf> vs; // 保存所有平面的数组
vector<lines> vls; // 保存所有初始反射路径（含碰撞点信息）的数组；
double moveLen = 1.0, width=2.0, swidth=0.4;
bool decrease = false;
double percent = 0.9;
GLfloat rValue = 1.0f, gValue = 1.0f, bValue = 1.0f;
double boxWidth = 1.44, boxLength = 1.44, boxHeight = 1.76;
CH3D hull; //全局变量，所有的空间凸包信息存储在这里；






/*---------------------我---是---分---割---线-----------（所有函数定义）-------------*/


/*~~~~~~~~~~~glut辅助函数~~~~~~~~~~~~~~~~*/

//初始化函数
void init(void)
{
	glEnable(GL_MAP2_VERTEX_3);
	glEnable(GL_DEPTH_TEST);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);//混合函数
	glEnable(GL_BLEND);
}

//捕捉鼠标移动函数
void motion(int x, int y)
{
	GLint dx,dy; //offset of mouse;

	dx = (GLint)x-mx;
	dy = (GLint)y-my;

	if(m_state == 0)
	{
		y_angle += dx*0.1f;
		x_angle += dy*0.1f;
	}
	else if(m_state == 1)
		dist += (dx+dy)*0.01f;
	
	mx = (GLint)x;
	my = (GLint)y;

	glutPostRedisplay();
}

//重绘制函数
void reshape(int w, int h)
{
	glViewport(0, 0, (GLint)w, (GLint)h);
}

//键盘设置函数
void keyboard(unsigned char key, int x, int y)
{
	switch(key)
	{
		case '0':
			m_state = 0;
			break;
		case '1':
			m_state = 1;
			break;
		default:
			break;
	}
}

//鼠标移动函数
void mouse(int button, int state, int x, int y)
{
	mx = (GLint)x;
	my = (GLint)y;
}




/*~~~~~~~~~~~功能函数~~~~~~~~~~~~~~~~*/

// 从所有出射光线里筛选从窗口入射的直线函数
void selectPoints() {
	for (int i=0; i<vp.size(); i++) {
		double k = vp[i].x / moveLen;
		if(vp[i].x > moveLen && vp[i].y/k>-swidth/2 && vp[i].y/k<swidth/2 && vp[i].z/k>-swidth/2 && vp[i].z/k<swidth/2) {
			Point p1(0,0,0);
			struct lines l = {};
			l.inner.push_back(p1);
			l.inner.push_back(vp[i]);
			vls.push_back(l);
		}
	}	
}

//空间两个点比较函数
bool cmp(Point a, Point b) {
	if(fabs(a.x - b.x) > threHold) return a.x < b.x;
	else if(fabs(a.y - b.y) > threHold) return a.y < b.y;
	else return a.z < b.z;
}

//绘制两点直线函数
void drawLines(Point p1, Point p2) {
	glColor3f( rValue, gValue, bValue);
	glLineWidth(lineWidth);
	glBegin(GL_LINES);
		glVertex3f( p1.x, p1.y, p1.z );
		glVertex3f(  p2.x, p2.y, p2.z );
	glEnd();
}

//绘制一个点函数
void drawPoint(Point p) {
	glPointSize(pointSize);
	glColor3f( rValue, gValue, bValue);
	glBegin(GL_POINTS);
		glVertex3f(p.x, p.y, p.z);
	glEnd();
}

// 计算空间s点到直线ab的距离；
double DistanceOfPointToLine(Point* a, Point* b, Point* s) 
{ 
	double ab = sqrt(pow((a->x - b->x), 2.0) + pow((a->y - b->y), 2.0) + pow((a->z - b->z), 2.0));
	double as = sqrt(pow((a->x - s->x), 2.0) + pow((a->y - s->y), 2.0) + pow((a->z - s->z), 2.0));
	double bs = sqrt(pow((s->x - b->x), 2.0) + pow((s->y - b->y), 2.0) + pow((s->z - b->z), 2.0));
	double cos_A = (pow(as, 2.0) + pow(ab, 2.0) - pow(bs, 2.0)) / (2 * ab*as);
	double sin_A = sqrt(1 - pow(cos_A, 2.0));
	return as*sin_A; 
}

// 计算空间两个点a，b的空间距离；
double distanceOfTwoPoints(Point a, Point b) {
	return pow((a.x-b.x), 2) + pow((a.z-b.z), 2) + pow((a.z-b.z), 2);
}

void drawTriangle(Point a, Point b, Point c) {
	glColor3f( 1.0f, 0.0f, 1.0f );
	//glColor4f( 0.5f, 0.5f, 0.0f,0.5f );
	glBegin(GL_TRIANGLES);
		glVertex3f(a.x, a.y, a.z);
		glVertex3f(b.x, b.y, b.z);
		glVertex3f(c.x, c.y, c.z);
	glEnd();
	drawLines(a, b);
	drawLines(b, c);
	drawLines(a, c);
}



/*~~~~~~~~~~~数据准备部分~~~~~~~~~~~~~~~~*/


//数据准备函数
void dataPrepare(){
	initSurfs(); //存入初始环境平面信息
	//initSurfs2(); // 换其他的平面场景
	//selectPoints(); //筛选入射直线的点信息
	initMultipleLight(); //产生初始所有出射光线数组的函数
	analyze(); // 得出了所有的反射经过点的信息；
	explicite(); //对在平面内的折射点进行包络体求解的算法
}

//测试入射线的函数
void initMultipleLight() {
	double k = 4; //最长可以辐射多远，比例参数
	Point p1(0,0.8,0.9);
	Point p2;
	int n = 20;
	double step = int(swidth*100/n)/100.0;
	for (double i=-swidth/2; i<=swidth/2; i+=step) {
		for (double j=-swidth/2; j<=swidth/2; j+=step) {
			Point p2((k+1)*moveLen-k*p1.x, (k+1)*i-k*p1.y, (k+1)*j-k*p1.z);
			struct lines l = {};
			l.inner.push_back(p1);
			l.inner.push_back(p2);
			vls.push_back(l);
		}
	}
}

//产生单根入射线的方法
void initSingleLight() {
	double k = 6;
	Point p1(0,0.3,0);
	Point p2((k+1)*moveLen-k*p1.x, (k+1)*-0.1-k*p1.y, (k+1)*0-k*p1.z);
	struct lines l = {};
	l.inner.push_back(p1);
	l.inner.push_back(p2);
	vls.push_back(l);
}


//初始化环境平面信息；
/*
boxWidth = 0.72, boxLength = 0.72, boxHeight = 0.88;
*/
void initSurfs() {
	struct surf sright = {boxWidth+moveLen,boxWidth+moveLen,-boxHeight/2,boxHeight/2,-boxLength/2,boxLength/2,1,0,0,-boxWidth-moveLen};
	vs.push_back(sright);
	struct surf sup = {moveLen,boxWidth+moveLen,boxHeight/2,boxHeight/2,-boxLength/2,boxLength/2,0,1,0,-boxHeight/2};
	vs.push_back(sup);
	struct surf sdown = {moveLen,boxWidth+moveLen,-boxHeight/2,-boxHeight/2,-boxLength/2,boxLength/2,0,1,0,boxHeight/2};
	vs.push_back(sdown);
	struct surf sfront = {moveLen,boxWidth+moveLen,-boxHeight/2,boxHeight/2,boxLength/2,boxLength/2,0,0,1,-boxLength/2};
	vs.push_back(sfront);
	struct surf sback = {moveLen,boxWidth+moveLen,-boxHeight/2,boxHeight/2,-boxLength/2,-boxLength/2,0,0,1,boxLength/2};
	vs.push_back(sback);
	struct surf sleft = {moveLen,moveLen,-boxHeight/2,boxHeight/2,-boxLength/2,boxLength/2,1,0,0,-moveLen, true,moveLen,moveLen,-swidth/2,swidth/2,-swidth/2,swidth/2};
	vs.push_back(sleft);
}


//循环折射函数
/*
让vls数组里每一条出射光线都在壳体中反复折射的主题函数
*/
void analyze() {
	int c = 3;
	for (int i=0; i<vls.size(); i++) {
		int index = 0;
		Point hp, next;
		struct lines line = vls[i];
		while (index < vls[i].inner.size()-1) {
			if (!collision(vls[i].inner[index], vls[i].inner[index+1], &hp, &next, i)) {
				vls[i].inner[index+1].x=hp.x;
				vls[i].inner[index+1].y=hp.y;
				vls[i].inner[index+1].z=hp.z;
				vls[i].inner.push_back(next);
			}
			index ++;
		}
	}
}

//对在平面内的折射点进行包络体求解的算法
void explicite() {
	for(int i=0; i<vs.size(); i++) {
		Point n1(vs[i].a, vs[i].b, vs[i].c);
		dealWith(vs[i].allPoints, vs[i].polygon, n1);
	}
	sortAndConstructThePolygonPoints();
}

//对包络平面的点进行进一步处理函数
void sortAndConstructThePolygonPoints() {
	vector<Point> allPolygonEdgePoints;
	for(int i=0; i<vs.size(); i++) {
		for(int j=0; j<vs[i].polygon[0].size(); j++)
			allPolygonEdgePoints.push_back(vs[i].polygon[0][j]);
		for(int j=1; j<vs[i].polygon[1].size()-1; j++)
			allPolygonEdgePoints.push_back(vs[i].polygon[1][j]);
		sort(vs[i].polygon[0].begin(), vs[i].polygon[0].end(), cmp);
		sort(vs[i].polygon[1].begin(), vs[i].polygon[1].end(), cmp);
	}
	hull.n = allPolygonEdgePoints.size();
	for(int i=0; i<hull.n; i++) {
		hull.P[i].x = allPolygonEdgePoints[i].x;
		hull.P[i].y = allPolygonEdgePoints[i].y;
		hull.P[i].z = allPolygonEdgePoints[i].z;
	}
	hull.create();
	polygons();
}

void polygons() {
	for(int i=0; i<hull.count.size(); i++) {
		if(hull.count[i].size() < 4) continue;
		else {
			vector<Point> allPoints;
			auto it=hull.count[i].begin();
			for(; it!=hull.count[i].end(); it++)
				allPoints.push_back(hull.P[*it]);
			it=hull.count[i].begin();
			Point a = hull.P[*it++];
			Point b = hull.P[*it++];
			Point c = hull.P[*it];
			Point n1 = (b-a)*(c-a);
			dealWith(allPoints, hull.polygons[i], n1);
		}
	}
}


/*~~~~~~~~~~~循环折射计算部分~~~~~~~~~~~~~~~~*/

/* 碰撞检验函数
*  有碰撞返回0，并把交叉点信息保存在hp里面，镜像点信息保存在next里，否则返回1；
*  p1是第一个点，p2是第二个点，hp是交叉点，next是下一个反射点，
*/
int collision(Point p1, Point p2, Point *hp, Point *next, int index) {
	for (int i=0; i<vs.size(); i++) {
		struct surf temp = vs[i];
		double result = (temp.a * p1.x + temp.b * p1.y + temp.c * p1.z + temp.d) * (temp.a * p2.x + temp.b * p2.y + temp.c * p2.z + temp.d);
		
		if (result <= 0) {
			int size = vls[index].surfIndex.size();
			if (result == 0) {
				if (size >=1 && vls[index].surfIndex[size-1] == i) continue;
				else if (size >=2 && vls[index].surfIndex[size-2] == i) continue;
			} 
			Point cross = getCross(p1, p2, temp);

			if(vls[index].surfIndex.size() > 0 && i == vls[index].surfIndex.back()) continue;
			
			if (validate(cross, temp)) {
					if (temp.isGate && cross.x >= temp.gxmin && cross.x <= temp.gxmax && 
						cross.y >= temp.gymin && cross.y <= temp.gymax &&
						cross.z >= temp.gzmin && cross.z <= temp.gzmax) continue;
				hp->x = cross.x;
				hp->y = cross.y;
				hp->z = cross.z;
				double k = -(temp.a*p2.x+temp.b*p2.y+temp.c*p2.z+temp.d)/(temp.a*temp.a+temp.b*temp.b+temp.c*temp.c);
				next->x=2*temp.a*k+p2.x;
				next->y=2*temp.b*k+p2.y;
				next->z=2*temp.c*k+p2.z;

				if(decrease) {
					next->x = percent * next->x + (1-percent)* hp->x;
					next->y = percent * next->y + (1-percent)* hp->y;
					next->z = percent * next->z + (1-percent)* hp->z;
				}
				vls[index].surfIndex.push_back(i);
				vs[i].allPoints.push_back(*hp);
				return 0;
			}			
		}
	}
	return 1;
}

//计算p1和p2与面域s之间的交点坐标函数，返回交点 temp
Point getCross(Point p1, Point p2, surf s) {
	Point temp;
	double k = (s.a*p1.x+s.b*p1.y+s.c*p1.z+s.d)/(s.a*(p1.x-p2.x)+s.b*(p1.y-p2.y)+s.c*(p1.z-p2.z));
	temp.x=p1.x+k*(p2.x-p1.x);
	temp.y=p1.y+k*(p2.y-p1.y);
	temp.z=p1.z+k*(p2.z-p1.z);
	return temp;
}

/**
* 验证交点是否在某平面之内的函数；
**/
int validate(Point cross, surf temp) {
	if(temp.xmin == temp.xmax && abs(cross.x - temp.xmin) > threHold) return 0;
	if(temp.ymin == temp.ymax && abs(cross.y - temp.ymin) > threHold) return 0;
	if(temp.zmin == temp.zmax && abs(cross.z - temp.zmin) > threHold) return 0;
	if(temp.xmin != temp.xmax && (cross.x < temp.xmin || cross.x > temp.xmax)) return 0;
	if(temp.ymin != temp.ymax && (cross.y < temp.ymin || cross.y > temp.ymax)) return 0;
	if(temp.zmin != temp.zmax && (cross.z < temp.zmin || cross.z > temp.zmax)) return 0;
	return 1;
}





/*~~~~~~~~~~~包络体计算部分~~~~~~~~~~~~~~~~*/

/*求平面内的最大包络多边形
参数解释：平面内所有点信息，用于存储多边形上下两半的二维数组，平面的法向量
*/
void dealWith(vector<Point> &allPoints, vector<Point> polygon[2], Point n1) {
	if(allPoints.size() < 2) return;
	Point a, b; //最小和最大两个极端顶点；
	a.x = allPoints[0].x;
	a.y = allPoints[0].y;
	a.z = allPoints[0].z;
	b.x = allPoints[0].x;
	b.y = allPoints[0].y;
	b.z = allPoints[0].z;
	for(int i=1; i<allPoints.size(); i++) {
		if(a.x - allPoints[i].x > eps) {
			a.x = allPoints[i].x;
			a.y = allPoints[i].y;
			a.z = allPoints[i].z;
		} else if(fabs(a.x - allPoints[i].x) < eps) {
			if(a.y - allPoints[i].y > eps) {
				a.x = allPoints[i].x;
				a.y = allPoints[i].y;
				a.z = allPoints[i].z;
			} else if(fabs(a.y - allPoints[i].y) < eps) {
				if(a.z - allPoints[i].z > eps) {
					a.x = allPoints[i].x;
					a.y = allPoints[i].y;
					a.z = allPoints[i].z;
				} 
			}
		}

		if(allPoints[i].x - b.x > eps) {
			b.x = allPoints[i].x;
			b.y = allPoints[i].y;
			b.z = allPoints[i].z;
		} else if(fabs(b.x - allPoints[i].x) < eps) {
			if(allPoints[i].y - b.y > eps) {
				b.x = allPoints[i].x;
				b.y = allPoints[i].y;
				b.z = allPoints[i].z;
			} else if(fabs(b.y - allPoints[i].y) < eps) {
				if(allPoints[i].z - b.z > eps) {
					b.x = allPoints[i].x;
					b.y = allPoints[i].y;
					b.z = allPoints[i].z;
				} 
			}
		}
	}

	if (fabs(a.x - b.x) + fabs(a.y - b.y) + fabs(a.z - b.z) < eps) {
		polygon[0].push_back(a);
		printf("两极值点相距过近，返回了直接");
		return;
	}

	polygon[0].push_back(a);
	polygon[0].push_back(b);
	polygon[1].push_back(a);
	polygon[1].push_back(b);
	vector<Point> p1, p2; // p1是直线左边所有点集合，p2是直线右边所有点集合
	Point mid ((a.x+b.x)/2, (a.y+b.y)/2, (a.z+b.z)/2); // 线段中点
	Point n2 (b.x-a.x, b.y-a.y, b.z-a.z); //两个极值点的线段所在的向量
	Point n3 (n1.y*n2.z-n2.y*n1.z, n2.x*n1.z-n1.x*n2.z, n1.x*n2.y-n2.x*n1.y); // 计算所在平面内的线段的法向量
	for (int i = 0; i < allPoints.size(); ++i)
	{
		Point temp (allPoints[i].x-mid.x, allPoints[i].y-mid.y, allPoints[i].z-mid.z); //点集合中任意一个点到直线中点的向量
		double value = n3.x*temp.x + n3.y*temp.y + n3.z*temp.z; //向量和平面内直线法向量的点积
		if(value > eps) p1.push_back(allPoints[i]);
		else if(value < -eps) p2.push_back(allPoints[i]);
	}
	FindPoint2(p1, a, b, mid, polygon[0], n1);
	FindPoint2(p2, a, b, mid, polygon[1], n1);
}

//平面求包主题算法
void FindPoint2(vector<Point> &p, Point a, Point b, Point mid, vector<Point> &polygon, Point &n) {
	if (p.size() == 0)
		return;
	Point pmax;
	pmax.x = p[0].x;
	pmax.y = p[0].y;
	pmax.z = p[0].z;
	double k, d;
	k = (b.y - a.y) / (b.x - a.x);
	d = a.y - k * a.x;
	double maxDis = DistanceOfPointToLine(&a, &b, &pmax), maxMid = distanceOfTwoPoints(pmax, mid);
	double newdist;
	for (int i = 1; i < p.size(); ++i)
	{
		newdist = DistanceOfPointToLine(&a, &b, &p[i]);
		if (newdist - maxDis > eps)
		{
			pmax.x = p[i].x;
			pmax.y = p[i].y;
			pmax.z = p[i].z;
			maxDis = newdist;
		}
		else if (fabs(newdist - maxDis) < eps)
		{	//选择距离线段ab中点最近的那个
			double dis1 = distanceOfTwoPoints(p[i], mid);
			if (dis1 < maxMid)
			{
				pmax.x = p[i].x;
				pmax.y = p[i].y;
				pmax.z = p[i].z;
				maxMid = dis1;
			}
		}
	}
	polygon.push_back(pmax);

	Point mid1 ((pmax.x+a.x)/2, (pmax.y+a.y)/2, (pmax.z+a.z)/2);
	Point mid2 ((pmax.x+b.x)/2, (pmax.y+b.y)/2, (pmax.z+b.z)/2);
	Point v1 (mid1.x-mid.x, mid1.y-mid.y, mid1.z-mid.z);
	Point v2 (mid2.x-mid.x, mid2.y-mid.y, mid2.z-mid.z);
	Point l1 (pmax.x-a.x, pmax.y-a.y, pmax.z-a.z); //两个极值点的线段所在的向量
	Point n1 (n.y*l1.z-l1.y*n.z, l1.x*n.z-n.x*l1.z, n.x*l1.y-l1.x*n.y); // 计算所在平面内的线段的法向量
	Point l2 (pmax.x-b.x, pmax.y-b.y, pmax.z-b.z); //两个极值点的线段所在的向量
	Point n2 (n.y*l2.z-l2.y*n.z, l2.x*n.z-n.x*l2.z, n.x*l2.y-l2.x*n.y); // 计算所在平面内的线段的法向量
	if(v1.x*n1.x+v1.y*n1.y+v1.z*n1.z < -eps) {
		n1.x *= -1;
		n1.y *= -1;
		n1.z *= -1;
	}
	double len = sqrt(n1.x*n1.x+n1.y*n1.y+n1.z*n1.z);
	n1.x /= len;
	n1.y /= len;
	n1.z /= len;
	if(v2.x*n2.x+v2.y*n2.y+v2.z*n2.z < -eps) {
		n2.x *= -1;
		n2.y *= -1;
		n2.z *= -1;
	}
	len = sqrt(n2.x*n2.x+n2.y*n2.y+n2.z*n2.z);
	n2.x /= len;
	n2.y /= len;
	n2.z /= len;

	/* 找出各自符合满足 Pmax,Pa 和 Pmax,Pb 的点 */
	vector<Point> p1, p2;
	for (int i = 0; i < p.size(); ++i)
	{
		Point temp1 (p[i].x-mid1.x, p[i].y-mid1.y, p[i].z-mid1.z);
		double value = temp1.x*n1.x+temp1.y*n1.y+temp1.z*n1.z;
		if(value > eps) p1.push_back(p[i]);
		else {
			Point temp2 (p[i].x-mid2.x, p[i].y-mid2.y, p[i].z-mid2.z);
			value = temp2.x*n2.x+temp2.y*n2.y+temp2.z*n2.z;
			if(value > eps) p2.push_back(p[i]);
		}

	}
	/* 递归寻找Pmax */
	FindPoint2(p1, pmax, a, mid1, polygon, n);
	FindPoint2(p2, pmax, b, mid2, polygon, n);
}







/*~~~~~~~~~~~主体逻辑函数~~~~~~~~~~~~~~~~*/

// 主函数
int main(int argc,char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(1200,800);
	glutInitWindowPosition(100,0);
	glutCreateWindow("RFID-model views");
	init();

	dataPrepare(); //准备一些数据

	printf("0 keydown means control the angle of the eye\n");
	printf("1 keydown means control the distance of the eye\n");

	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutMainLoop();	
	return 0;
}

//主显示函数
void display(void)
{
	int i, j;
	GLint rect[4];
	glClearColor(0.0f,0.0f,0.0f,0.0f);
	glClearDepth(1.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glGetIntegerv(GL_VIEWPORT, rect);
	if(rect[3] < 1) rect[3]=1;
	gluPerspective(30.0, 1.0*rect[2]/rect[3], 0.1, 20.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslated(0.0, 0.0, -dist);
	glRotatef(x_angle, 1.0f, 0.0f, 0.0f);
	glRotatef(y_angle, 0.0f, 1.0f, 0.0f);

	//显示坐标轴及周围环境
	drawRefections();	
	drawEnvironments();	
	//drawRefections();	
	glutSwapBuffers();
}





/*~~~~~~~~~~~环境情况绘制~~~~~~~~~~~~~~~~*/

//绘制环境信息
void drawEnvironments() {
	drawCoordinates(); // 绘制坐标轴
	drawSurfaces(); // 绘制平面信息
	drawEnvironmentLines(); // 绘制平面的轮廓直线
}


//绘制坐标轴函数
void drawCoordinates(void)
{
	glColor3f(1.0f,0.0f,0.0f); //画红色的x轴
	glBegin(GL_LINES);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(1.0f, 0.0f, 0.0f);
	glEnd();
	glColor3f(0.0,1.0,0.0); //画绿色的y轴
	glBegin(GL_LINES);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 1.0f, 0.0f);
	glEnd();
	glColor3f(0.0,0.0,1.0); //画蓝色的z轴
	glBegin(GL_LINES);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 1.0f);
	glEnd();

}

// 绘制环境平面信息
/*
boxWidth = 0.72, boxLength = 0.72, boxHeight = 0.88;
*/
void drawSurfaces() {
	
	glColor4f( 0.5f, 0.5f, 0.5f,0.5f );
	//绘制背后面
	glBegin( GL_QUADS );
		glVertex3f( boxWidth+moveLen, boxHeight/2, -boxLength/2 );
		glVertex3f( boxWidth+moveLen, boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, -boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, -boxHeight/2, -boxLength/2 );
	glEnd();
	glColor4f( 0.8f, 0.8f, 0.8f,0.5f );
	//绘制上面
	glBegin( GL_QUADS );
		glVertex3f( moveLen, boxHeight/2, -boxLength/2 );
		glVertex3f( moveLen, boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, boxHeight/2, -boxLength/2 );
	glEnd();
	//绘制下面
	glBegin( GL_QUADS );
		glVertex3f( moveLen, -boxHeight/2, -boxLength/2 );
		glVertex3f( moveLen, -boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, -boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, -boxHeight/2, -boxLength/2 );
	glEnd();
	////绘制前面
	glBegin( GL_QUADS );
		glVertex3f( moveLen, boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, -boxHeight/2, boxLength/2 );
		glVertex3f( moveLen, -boxHeight/2, boxLength/2 );
	glEnd();
	////绘制背后面
	glBegin( GL_QUADS );
		glVertex3f( moveLen, boxHeight/2, -boxLength/2 );
		glVertex3f( boxWidth+moveLen, boxHeight/2, -boxLength/2 );
		glVertex3f( boxWidth+moveLen, -boxHeight/2, -boxLength/2 );
		glVertex3f( moveLen, -boxHeight/2, -boxLength/2 );
	glEnd();

	//绘制正面
	glColor4f( 1.0f,1.0f,1.0f,0.5f);
	glBegin( GL_QUADS );
		glVertex3f( moveLen, -boxHeight/2, -boxLength/2 );
		glVertex3f( moveLen, -swidth/2, -boxLength/2 );
		glVertex3f( moveLen, -swidth/2, boxLength/2 );
		glVertex3f( moveLen, -boxHeight/2, boxLength/2 );
	glEnd();

	glBegin( GL_QUADS );
		glVertex3f( moveLen, boxHeight/2, -boxLength/2 );
		glVertex3f( moveLen, swidth/2, -boxLength/2 );
		glVertex3f( moveLen, swidth/2, boxLength/2 );
		glVertex3f( moveLen, boxHeight/2, boxLength/2 );
	glEnd();

	glBegin( GL_QUADS );
		glVertex3f( moveLen, -swidth/2 , boxLength/2);
		glVertex3f( moveLen, -swidth/2 , swidth/2);
		glVertex3f( moveLen, swidth/2 , swidth/2);
		glVertex3f( moveLen, swidth/2 , boxLength/2);
	glEnd();

	glBegin( GL_QUADS );
		glVertex3f( moveLen, -swidth/2 , -boxLength/2);
		glVertex3f( moveLen, -swidth/2 , -swidth/2);
		glVertex3f( moveLen, swidth/2 , -swidth/2);
		glVertex3f( moveLen, swidth/2 , -boxLength/2);
	glEnd();

}

//绘制立方体边框直线
/*
x:boxWidth = 0.72, z:boxLength = 0.72, y:boxHeight = 0.88;
*/
void drawEnvironmentLines(){
	glColor3f( 0.0, 0.0, 0.0 );
	//lineWidth = 1.0;
	glLineWidth(lineWidth);
	glBegin(GL_LINES);
		glVertex3f( moveLen, -boxHeight/2, -boxLength/2 );
		glVertex3f( moveLen, -boxHeight/2, boxLength/2 );
		glVertex3f( moveLen, -boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, -boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, -boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, -boxHeight/2, -boxLength/2 );
		glVertex3f( boxWidth+moveLen, -boxHeight/2, -boxLength/2 );
		glVertex3f( moveLen, -boxHeight/2, -boxLength/2 );

		glVertex3f( moveLen, boxHeight/2, -boxLength/2 );
		glVertex3f( moveLen, boxHeight/2, boxLength/2 );
		glVertex3f( moveLen, boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, boxHeight/2, -boxLength/2 );
		glVertex3f( boxWidth+moveLen, boxHeight/2, -boxLength/2 );
		glVertex3f( moveLen, boxHeight/2, -boxLength/2 );

		glVertex3f( moveLen, -boxHeight/2, -boxLength/2 );
		glVertex3f( moveLen, boxHeight/2, -boxLength/2 );
		glVertex3f( moveLen, -boxHeight/2, boxLength/2 );
		glVertex3f( moveLen, boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, -boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, boxHeight/2, boxLength/2 );
		glVertex3f( boxWidth+moveLen, -boxHeight/2, -boxLength/2 );
		glVertex3f( boxWidth+moveLen, boxHeight/2, -boxLength/2 );


		//内框
		glVertex3f( moveLen, -swidth/2 , -swidth/2);
		glVertex3f( moveLen, -swidth/2 , swidth/2);
		glVertex3f( moveLen, -swidth/2 , swidth/2);
		glVertex3f( moveLen, swidth/2 , swidth/2);
		glVertex3f( moveLen, swidth/2 , swidth/2);
		glVertex3f( moveLen, swidth/2 , -swidth/2);
		glVertex3f( moveLen, swidth/2 , -swidth/2);
		glVertex3f( moveLen, -swidth/2 , -swidth/2);

	glEnd();
}





/*~~~~~~~~~~~折射情况绘制~~~~~~~~~~~~~~~~*/


//后期绘制折射所有信息
void drawRefections() {
	//drawRefectionTracks(); //绘制所有入射光线的折射路径
	//drawRefectionPoints(); //绘制所有光线的所有折射点
	//drawReflectionPointsInSurfs(); //绘制所有折射路径和所有平面的有效交点
	//drawRefectionPolygonLines(); //仅绘制在平面上的包络体的最外形轮廓
	//drawPolygonSurfs();
	drawSpacePolygonFaces();
}


//绘制所有入射光线的折射路径
void drawRefectionTracks() {
	rValue = 1.0f;
	gValue = 1.0f;
	bValue = 0.0f;
	lineWidth = 0.1f;
	for (int i=0; i<vls.size(); i++) {
		for(int j=0; j<vls[i].inner.size()-1; j++) {
			drawLines(vls[i].inner[j], vls[i].inner[j+1]);
		}
	}
}

//绘制所有光线的所有关键折射点
void drawRefectionPoints() {
	rValue = 1.0f;
	gValue = 1.0f;
	bValue = 0.0f;
	pointSize = 2.0f;
	for (int i=0; i<vls.size(); i++) {
		for(int j=0; j<vls[i].inner.size(); j++) {
			drawPoint(vls[i].inner[j]);
		}
	}
}

//仅绘制在平面上的包络体的最外形轮廓
void drawRefectionPolygonLines() {
	rValue = 1.0f;
	gValue = 1.0f;
	bValue = 0.0f;
	lineWidth = 1.0f;
	for(int i=0; i<vs.size(); i++) {
		sort(vs[i].polygon[0].begin(), vs[i].polygon[0].end(), cmp);
		sort(vs[i].polygon[1].begin(), vs[i].polygon[1].end(), cmp);
		for(int j=0; j<vs[i].polygon[0].size()-1; j++) {
			drawLines(vs[i].polygon[0][j], vs[i].polygon[0][j+1]);
		}
		for(int j=0; j<vs[i].polygon[1].size()-1; j++) {
			drawLines(vs[i].polygon[1][j], vs[i].polygon[1][j+1]);
		}
	}
}

//绘制所有折射路径和所有平面的有效交点
void drawReflectionPointsInSurfs() {
	rValue = 1.0f;
	gValue = 1.0f;
	bValue = 0.0f;
	pointSize = 2.0f;
	for (int i=0; i<vs.size(); i++) {
		for(int j=0; j<vs[i].allPoints.size(); j++) {
			drawPoint(vs[i].allPoints[j]);
		}
	}
}

//绘制在内壳上的最大包络平面（用红色高亮）
void drawPolygonSurfs() {
	for(int i=0; i<vs.size(); i++) {
		glBegin(GL_POLYGON);
			glColor3f(1.0f,0.0f,0.0f);
			for(int j=0; j<vs[i].polygon[0].size(); j++)
				glVertex3f(vs[i].polygon[0][j].x, vs[i].polygon[0][j].y, vs[i].polygon[0][j].z);
			for(int j=0; j<vs[i].polygon[1].size(); j++)
				glVertex3f(vs[i].polygon[1][j].x, vs[i].polygon[1][j].y, vs[i].polygon[1][j].z);
		glEnd();
	}

}

void drawSpacePolygonFaces() {
	for(int i=0; i<hull.count.size(); i++) {
		if (hull.count[i].size() < 4) {
			auto it=hull.count[i].begin();
			Point a = hull.P[*it++];
			Point b = hull.P[*it++];
			Point c = hull.P[*it];
			drawTriangle(a, b, c);
		} else {
			sort(hull.polygons[i][0].begin(), hull.polygons[i][0].end(), cmp);
			sort(hull.polygons[i][1].begin(), hull.polygons[i][1].end(), cmp);
			for(int j=0; j<hull.polygons[i][0].size()-1; j++) {
				drawLines(hull.polygons[i][0][j], hull.polygons[i][0][j+1]);
			}
			for(int j=0; j<hull.polygons[i][1].size()-1; j++) {
				drawLines(hull.polygons[i][1][j], hull.polygons[i][1][j+1]);
			}
			glBegin(GL_POLYGON);
				glColor3f(1.0f,0.0f,1.0f);
				for(int j=0; j<hull.polygons[i][0].size(); j++)
					glVertex3f(hull.polygons[i][0][j].x, hull.polygons[i][0][j].y, hull.polygons[i][0][j].z);
				for(int j=0; j<hull.polygons[i][1].size(); j++)
					glVertex3f(hull.polygons[i][1][j].x, hull.polygons[i][1][j].y, hull.polygons[i][1][j].z);
			glEnd();
		}
	}
}
