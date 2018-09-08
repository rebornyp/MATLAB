#include <stdio.h>
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <glut.h>
#include <vector>
#include <algorithm>
using namespace std;

#define pi 3.1415926
#define threHold 0.001

typedef struct point {
	double x, y, z;
}Point;

struct surf {
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double a,b,c,d;
	bool isGate;
	double gxmin, gxmax, gymin, gymax, gzmin, gzmax;
	vector<point> allPoints;
	vector<point> polygon;
};

struct lines {
	vector<point> inner;
	vector<int> surfIndex;
};

void init(void);
void reshape(int w, int h);
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void display(void);
void drawCoordinates(void);
void drawTetrahedron(void);
void drawLines();
void dataPrepare();
void drawSurfaces();
void testPoints();
void testSurfs(int);
void drawTest();
void drawRefections();
void initMultipleLight();
void initPoints();
void analyze();
int collision(point p1, point p2, struct point *hp, struct point *next, int index);
point getCross(point p1, point p2, surf s);
void drawPoint(point p);
void drawLines(point p1, point p2);
void drawEnvironments();
void drawSurfaces();
void initSurfs();
void createSphereArray(int num, double r);
void initSingleLight();
int validate(point cross, surf temp);
void explicite();
void dealWith(vector<point> &allPoints, vector<point> &polygon, point n1);
void drawPolygons();
double DistanceOfPointToLine(Point* a, Point* b, Point* s);
double distanceOfTwoPoints(Point a, Point b);
void FindPoint2(vector<point> &p, Point a, Point b, Point mid, vector<point> &polygon);
void testSurfPolygons(vector<point> &p, vector<point> &polygon, point n1);
void testSurfPolygons2();
void drawRefectionTracks();
void drawRefectionPoints();
void drawRefectionPolygons();


GLint mx,my; //position of mouse
GLint m_state=0; //mouse usage
GLfloat x_angle=0.0f, y_angle=0.0f; //angle of eye
GLfloat dist=10.0f; //distance from the eye
double lineWidth = 0.01, pointSize = 0.1;
vector<point> vp, vtemp;
vector<surf> vs; // 保存所有平面的数组
vector<lines> vls; // 保存所有反射路径（含碰撞点信息）的数组；
double moveLen = 1.0, width=2.0, swidth=0.5;
bool decrease = false;
double percent = 0.9;
GLfloat rValue = 1.0f, gValue = 1.0f, bValue = 1.0f;


bool cmp(point a, point b) {
	if(a.x != b.x) return a.x < b.x;
	else return a.y < b.y;
}

// 主函数
int main(int argc,char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(1200,800);
	glutInitWindowPosition(100,0);
	glutCreateWindow("RFID-model views");
	init();
	//testSurfPolygons2();
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
	//drawTest();
	glutSwapBuffers();
}

//创立一开始球面的所有坐标个数及球半径信息
void createSphereArray(int num, double r) {
	double x,y,z, l1,l2,step1, step2;
	l2=0;
	l1=-90;
	step2=360/num;
	step1=180/num;

	for (int i=0; i<2*num; l2+=step2/2, i++) {
		x=r*cos(pi*l2/180);
		y=r*sin(pi*l2/180);
		for (int j=0; j<num; l1+=step1, j++) {
			double tx=x*cos(pi*l1/180);
			double ty=y*cos(pi*l1/180);
			z = r*sin(pi*l1/180);
			struct point p = {tx, ty, z};
			vp.push_back(p);
		}
	}	
}

//初始化函数
void init(void)
{
	glEnable(GL_MAP2_VERTEX_3);
	glEnable(GL_DEPTH_TEST);
	//glColor4f(1.0f, 1.0f, 1.0f, 0.5f);//颜色0.5 alpha值
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

//绘制两点直线函数
void drawLines(point p1, point p2) {
	glColor3f( rValue, gValue, bValue);
	glLineWidth(lineWidth);
	glBegin(GL_LINES);
		glVertex3f( p1.x, p1.y, p1.z );
		glVertex3f(  p2.x, p2.y, p2.z );
	glEnd();
}

//绘制一个点函数
void drawPoint(point p) {
	glPointSize(pointSize);
	glColor3f( rValue, gValue, bValue);
	glBegin(GL_POINTS);
		glVertex3f(p.x, p.y, p.z);
	glEnd();
}

//初始化环境平面信息；
void initSurfs() {
	struct surf sright = {width+moveLen,width+moveLen,-width/2,width/2,-width/2,width/2,1,0,0,-width-moveLen};
	vs.push_back(sright);
	struct surf sup = {moveLen,width+moveLen,width/2,width/2,-width/2,width/2,0,1,0,-width/2};
	vs.push_back(sup);
	struct surf sdown = {moveLen,width+moveLen,-width/2,-width/2,-width/2,width/2,0,1,0,width/2};
	vs.push_back(sdown);
	struct surf sfront = {moveLen,width+moveLen,-width/2,width/2,width/2,width/2,0,0,1,-width/2};
	vs.push_back(sfront);
	struct surf sback = {moveLen,width+moveLen,-width/2,width/2,-width/2,-width/2,0,0,1,width/2};
	vs.push_back(sback);
	struct surf sleft = {moveLen,moveLen,-width/2,width/2,-width/2,width/2,1,0,0,-moveLen, true,moveLen,moveLen,-swidth/2,swidth/2,-swidth/2,swidth/2};
	vs.push_back(sleft);

	// 三角片曲面信息；
	//struct surf triangle = {};

}

void initSurfs2() {

}

// 绘制环境平面信息
void drawSurfaces() {
	
	glColor4f( 0.5f, 0.5f, 0.5f,0.5f );
	//绘制背后面
	glBegin( GL_QUADS );
		glVertex3f( width+moveLen, width/2, -width/2 );
		glVertex3f( width+moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, -width/2 );
	glEnd();
	glColor4f( 0.8f, 0.8f, 0.8f,0.5f );
	//绘制上面
	glBegin( GL_QUADS );
		glVertex3f( moveLen, width/2, -width/2 );
		glVertex3f( moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, width/2, -width/2 );
	glEnd();
	//绘制下面
	glBegin( GL_QUADS );
		glVertex3f( moveLen, -width/2, -width/2 );
		glVertex3f( moveLen, -width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, -width/2 );
	glEnd();
	////绘制前面
	glBegin( GL_QUADS );
		glVertex3f( moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, width/2 );
		glVertex3f( moveLen, -width/2, width/2 );
	glEnd();
	////绘制背后面
	glBegin( GL_QUADS );
		glVertex3f( moveLen, width/2, -width/2 );
		glVertex3f( width+moveLen, width/2, -width/2 );
		glVertex3f( width+moveLen, -width/2, -width/2 );
		glVertex3f( moveLen, -width/2, -width/2 );
	glEnd();

	//绘制正面
	glColor4f( 1.0f,1.0f,1.0f,0.5f);
	glBegin( GL_QUADS );
		glVertex3f( moveLen, -width/2, -width/2 );
		glVertex3f( moveLen, -swidth/2, -width/2 );
		glVertex3f( moveLen, -swidth/2, width/2 );
		glVertex3f( moveLen, -width/2, width/2 );
	glEnd();

	glBegin( GL_QUADS );
		glVertex3f( moveLen, width/2, -width/2 );
		glVertex3f( moveLen, swidth/2, -width/2 );
		glVertex3f( moveLen, swidth/2, width/2 );
		glVertex3f( moveLen, width/2, width/2 );
	glEnd();

	glBegin( GL_QUADS );
		glVertex3f( moveLen, -swidth/2 , width/2);
		glVertex3f( moveLen, -swidth/2 , swidth/2);
		glVertex3f( moveLen, swidth/2 , swidth/2);
		glVertex3f( moveLen, swidth/2 , width/2);
	glEnd();

	glBegin( GL_QUADS );
		glVertex3f( moveLen, -swidth/2 , -width/2);
		glVertex3f( moveLen, -swidth/2 , -swidth/2);
		glVertex3f( moveLen, swidth/2 , -swidth/2);
		glVertex3f( moveLen, swidth/2 , -width/2);
	glEnd();

}

//绘制立方体边框直线
void drawLines(){
	glColor3f( 0.0, 0.0, 0.0 );
	glLineWidth(lineWidth);
	glBegin(GL_LINES);
		glVertex3f( moveLen, -width/2, -width/2 );
		glVertex3f( moveLen, -width/2, width/2 );
		glVertex3f( moveLen, -width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, -width/2 );
		glVertex3f( width+moveLen, -width/2, -width/2 );
		glVertex3f( moveLen, -width/2, -width/2 );

		glVertex3f( moveLen, width/2, -width/2 );
		glVertex3f( moveLen, width/2, width/2 );
		glVertex3f( moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, width/2, -width/2 );
		glVertex3f( width+moveLen, width/2, -width/2 );
		glVertex3f( moveLen, width/2, -width/2 );

		glVertex3f( moveLen, -width/2, -width/2 );
		glVertex3f( moveLen, width/2, -width/2 );
		glVertex3f( moveLen, -width/2, width/2 );
		glVertex3f( moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, width/2 );
		glVertex3f( width+moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, -width/2 );
		glVertex3f( width+moveLen, width/2, -width/2 );


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

//绘制环境信息
void drawEnvironments() {
	drawCoordinates(); // 绘制坐标轴
	drawSurfaces(); // 绘制平面信息
	drawLines(); // 绘制平面的轮廓直线

	drawPolygons();
}

void drawPolygons() {

}


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

void drawRefectionPolygons() {
	rValue = 1.0f;
	gValue = 0.0f;
	bValue = 0.0f;
	lineWidth = 5.0f;
	for(int i=0; i<vs.size(); i++) {
		for(int j=0; j<vs[i].polygon.size(); j++) {
			for(int k=j+1; k<vs[i].polygon.size(); k++)
				drawLines(vs[i].polygon[j], vs[i].polygon[k]);
		}
	}
}

//后期绘制折射线的函数
void drawRefections() {
	//drawRefectionTracks();
	drawRefectionPoints();
	drawRefectionPolygons();
}

// 从入口射入多跟入射光线；
void initPoints() {
	for (int i=0; i<vp.size(); i++) {
		double k = vp[i].x / moveLen;
		if(vp[i].x > moveLen && vp[i].y/k>-swidth/2 && vp[i].y/k<swidth/2 && vp[i].z/k>-swidth/2 && vp[i].z/k<swidth/2) {
			struct point p1={0,0,0};
			struct lines l = {};
			l.inner.push_back(p1);
			l.inner.push_back(vp[i]);
			vls.push_back(l);
		}
	}	
}

//数据准备函数
void dataPrepare(){
	//createSphereArray(20, 6); //绘制球面的所有点坐标
	initSurfs(); //存入初始环境平面信息
	//initSurfs2(); // 换其他的平面场景
	//initPoints(); //节选入射直线的点信息
	initMultipleLight();
	//initSingleLight();
	//testSurfs();
	//testPoints();
	analyze(); // 得出了所有的反射经过点的信息；
	explicite();

}

//测试入射线的函数
void initMultipleLight() {
	double k = 5;
	struct point p1={0,0.8,0.9};
	int n = 20;
	double step = int(0.5*100/n)/100.0;
	for (double i=-0.25; i<=0.25; i+=step) {
		for (double j=-0.25; j<=0.25; j+=step) {
			double px = (k+1)*moveLen-k*p1.x;
			double py = (k+1)*i-k*p1.y;
			double pz = (k+1)*j-k*p1.z;
			struct point p2={px,py,pz};
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
	struct point p1={0,0.3,0};
	struct point p2={(k+1)*moveLen-k*p1.x, (k+1)*-0.1-k*p1.y, (k+1)*0-k*p1.z};
	struct lines l = {};
	l.inner.push_back(p1);
	l.inner.push_back(p2);
	vls.push_back(l);
}

/**
* 验证交点是否在某平面之内的函数；
**/
int validate(point cross, surf temp) {
	//printf("和环境碰撞交点坐标是：(%f,%f,%f)\n", cross.x,cross.y,cross.z);
	//printf("碰撞面信息x：[%.2f %.2f]，y：[%.2f %.2f]，z：[%.2f %.2f]\n", temp.xmin,temp.xmax,temp.ymin,temp.ymax,temp.zmin,temp.zmax);
	//printf("welcome to come here....\n");
	if(temp.xmin == temp.xmax && abs(cross.x - temp.xmin) > threHold) return 0;
	//printf("1");
	if(temp.ymin == temp.ymax && abs(cross.y - temp.ymin) > threHold) return 0;
	//printf("2");
	if(temp.zmin == temp.zmax && abs(cross.z - temp.zmin) > threHold) return 0;
	//printf("3");
	if(temp.xmin != temp.xmax && (cross.x < temp.xmin || cross.x > temp.xmax)) return 0;
	//printf("4");
	if(temp.ymin != temp.ymax && (cross.y < temp.ymin || cross.y > temp.ymax)) return 0;
	//printf("5");
	if(temp.zmin != temp.zmax && (cross.z < temp.zmin || cross.z > temp.zmax)) return 0;
	//printf("6");
	return 1;
}

//计算p1和p2与面域s之间的交点坐标函数，返回交点 temp
point getCross(point p1, point p2, surf s) {
	struct point temp = {};
	double k = (s.a*p1.x+s.b*p1.y+s.c*p1.z+s.d)/(s.a*(p1.x-p2.x)+s.b*(p1.y-p2.y)+s.c*(p1.z-p2.z));
	temp.x=p1.x+k*(p2.x-p1.x);
	temp.y=p1.y+k*(p2.y-p1.y);
	temp.z=p1.z+k*(p2.z-p1.z);
	//temp.x = int(temp.x*100+0.5)/100.0;
	//temp.y = int(temp.y*100+0.5)/100.0;
	//temp.z = int(temp.z*100+0.5)/100.0;
	return temp;
}

/*
*  有碰撞返回0，并把交叉点信息保存在hp里面，镜像点信息保存在next里，否则返回1；
*  p1是第一个点，p2是第二个点，hp是交叉点，next是下一个反射点，
*/
int collision(point p1, point p2, struct point *hp, struct point *next, int index) {
	for (int i=0; i<vs.size(); i++) {
		struct surf temp = vs[i];
		double result = (temp.a * p1.x + temp.b * p1.y + temp.c * p1.z + temp.d) * (temp.a * p2.x + temp.b * p2.y + temp.c * p2.z + temp.d);
		
		if (result <= 0) {
			//printf("碰撞了。。。\n");
			int size = vls[index].surfIndex.size();
			if (result == 0) {
				if (size >=1 && vls[index].surfIndex[size-1] == i) continue;
				else if (size >=2 && vls[index].surfIndex[size-2] == i) continue;
			} 
			struct point cross = getCross(p1, p2, temp);

			//printf("和环境碰撞交点坐标是：(%.2f,%.2f,%.2f)\n", cross.x,cross.y,cross.z);
			//printf("和环境碰撞交点坐标是：(%f,%f,%f)\n", cross.x,cross.y,cross.z);
			
			//printf("碰撞面信息x：[%.2f %.2f]，y：[%.2f %.2f]，z：[%.2f %.2f]\n", temp.xmin,temp.xmax,temp.ymin,temp.ymax,temp.zmin,temp.zmax);
			
			//若是前后都同一平面则不允许
			if(vls[index].surfIndex.size() > 0 && i == vls[index].surfIndex.back()) continue;
			
			if (validate(cross, temp)) {
					//printf("碰撞交点是在碰撞面的总范围内了,并且碰撞面是");
					//printf(temp.isGate ? "出口曲面\n":"普通曲面\n");
					if (temp.isGate && cross.x >= temp.gxmin && cross.x <= temp.gxmax && 
						cross.y >= temp.gymin && cross.y <= temp.gymax &&
						cross.z >= temp.gzmin && cross.z <= temp.gzmax) continue;
				//printf("碰撞交点是在碰撞面的有效范围内了\n");
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

				//next->x=int(next->x*100+0.5)/100.0;
				//next->y=int(next->y*100+0.5)/100.0;
				//next->z=int(next->z*100+0.5)/100.0;
				//printf("镜像点是（%f %f %f）,碰撞面是第%d个平面\n", next->x, next->y,next->z,i);
				vls[index].surfIndex.push_back(i);
				vs[i].allPoints.push_back(*hp);
				return 0;
			}			
		}
	}
	return 1;
}

//循环折射函数
void analyze() {
	int c = 3;
	for (int i=0; i<vls.size(); i++) {
		//printf("第%d条直线...\n", i);
		int index = 0;
		struct point hp, next;
		struct lines line = vls[i];
		while (index < vls[i].inner.size()-1) {
			//if (index > c) break;
			//printf("第几个点：%d\n", index);
			/*printf("直线起始和终止两点分别是（%.2f, %.2f, %.2f），（%.2f, %.2f, %.2f）\n", vls[i].inner[index].x, vls[i].inner[index].y, vls[i].inner[index].z,
					vls[i].inner[index+1].x, vls[i].inner[index+1].y, vls[i].inner[index+1].z);*/
			if (!collision(vls[i].inner[index], vls[i].inner[index+1], &hp, &next, i)) {
				//printf("第%d个交叉点是（%f, %f, %f）\n", index, hp.x, hp.y, hp.z);
				vls[i].inner[index+1].x=hp.x;
				vls[i].inner[index+1].y=hp.y;
				vls[i].inner[index+1].z=hp.z;
				vls[i].inner.push_back(next);
				//vtemp.push_back(next);
			}
			index ++;
		}
	}

	/*for (int j=0; j<vls[0].surfIndex.size(); j++) {
		printf("%d\n", vls[0].surfIndex[j]);
	}*/
}

void explicite() {
	for(int i=0; i<vs.size(); i++) {
		point n1 = {vs[i].a, vs[i].b, vs[i].c};
		dealWith(vs[i].allPoints, vs[i].polygon, n1);
		testSurfs(i);
		//testSurfPolygons(vs[i].allPoints, vs[i].polygon, n1);
	}
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


void FindPoint2(vector<point> &p, Point a, Point b, Point mid, vector<point> &polygon) {
	if (p.size() == 0)
		return;
	/* 找出 Pmax 点 */
	Point pmax;
	pmax.x = p[0].x;
	pmax.y = p[0].y;
	pmax.z = p[0].z;
	double k, d;
	k = (b.y - a.y) / (b.x - a.x);
	d = a.y - k * a.x;
	double dist = DistanceOfPointToLine(&a, &b, &pmax);
	double newdist, maxDis = -1;
	for (int i = 1; i < p.size(); ++i)
	{
		newdist = DistanceOfPointToLine(&a, &b, &p[i]);
		if (newdist - dist > threHold)
		{
			pmax.x = p[i].x;
			pmax.y = p[i].y;
			pmax.z = p[i].z;
			maxDis = distanceOfTwoPoints(pmax, mid);
		}
		else if (fabs(newdist - dist) < threHold)
		{	//选择距离线段ab中点最近的那个
			double dis1 = distanceOfTwoPoints(p[i], mid);
			if (maxDis != -1 && dis1 < maxDis)
			{
				pmax.x = p[i].x;
				pmax.y = p[i].y;
				pmax.z = p[i].z;
			}
		}
	}
	printf("Pmax(%f, %f, %f)\n", pmax.x, pmax.y, pmax.z);
	polygon.push_back(pmax);

	Point mid1 = {(pmax.x+a.x)/2, (pmax.y+a.y)/2, (pmax.z+a.z)/2};
	Point mid2 = {(pmax.x+b.x)/2, (pmax.y+b.y)/2, (pmax.z+b.z)/2};
	Point v1 = {mid1.x-mid.x, mid1.y-mid.y, mid1.z-mid.z};
	Point v2 = {mid2.x-mid.x, mid2.y-mid.y, mid2.z-mid.z};

	/* 找出各自符合满足 Pmax,Pa 和 Pmax,Pb 的点 */
	vector<point> p1, p2;
	for (int i = 0; i < p.size(); ++i)
	{
		Point temp1 = {p[i].x-mid1.x, p[i].y-mid1.y, p[i].z-mid1.z};
		if(temp1.x*v1.x+temp1.y*v1.y+temp1.z*v1.z > threHold) p1.push_back(p[i]);
		else {
			Point temp2 = {p[i].x-mid2.x, p[i].y-mid2.y, p[i].z-mid2.z};
			if(temp2.x*v2.x+temp2.y*v2.y+temp2.z*v2.z > threHold) p2.push_back(p[i]);
		}

	}
	printf("v1.size=%d, v2.size=%d\n", p1.size(), p2.size());
	/* 递归寻找Pmax */
	FindPoint2(p1, pmax, a, mid1, polygon);
	printf("左边的算完了");
	FindPoint2(p2, pmax, b, mid2, polygon);
	printf("右边的算完了");
	//free(&p1);
	//free(&p2);
}


void dealWith(vector<point> &allPoints, vector<point> &polygon, point n1) {
	if(allPoints.size() < 2) return;
	Point a, b; //最小和最大两个极端顶点；
	a.x = allPoints[0].x;
	a.y = allPoints[0].y;
	a.z = allPoints[0].z;
	b.x = allPoints[0].x;
	b.y = allPoints[0].y;
	b.z = allPoints[0].z;
	for(int i=1; i<allPoints.size(); i++) {
		if(a.x - allPoints[i].x > threHold) {
			a.x = allPoints[i].x;
			a.y = allPoints[i].y;
			a.z = allPoints[i].z;
		} else if(fabs(a.x - allPoints[i].x) < threHold) {
			if(a.y - allPoints[i].y > threHold) {
				a.x = allPoints[i].x;
				a.y = allPoints[i].y;
				a.z = allPoints[i].z;
			} else if(fabs(a.y - allPoints[i].y) < threHold) {
				if(a.z - allPoints[i].z > threHold) {
					a.x = allPoints[i].x;
					a.y = allPoints[i].y;
					a.z = allPoints[i].z;
				} 
			}
		}

		if(allPoints[i].x - b.x > threHold) {
			b.x = allPoints[i].x;
			b.y = allPoints[i].y;
			b.z = allPoints[i].z;
		} else if(fabs(b.x - allPoints[i].x) < threHold) {
			if(allPoints[i].y - b.y > threHold) {
				b.x = allPoints[i].x;
				b.y = allPoints[i].y;
				b.z = allPoints[i].z;
			} else if(fabs(b.y - allPoints[i].y) < threHold) {
				if(allPoints[i].z - b.z > threHold) {
					b.x = allPoints[i].x;
					b.y = allPoints[i].y;
					b.z = allPoints[i].z;
				} 
			}
		}
	}

	if (fabs(a.x - b.x) + fabs(a.y - b.y) + fabs(a.z - b.z) < threHold) {
		polygon.push_back(a);
		printf("两极值点相距过近，返回了直接");
		return;
	}

	polygon.push_back(a);
	polygon.push_back(b);
	printf("Pmin(%f, %f, %f)\n", a.x, a.y, a.z);
	printf("Pmax(%f, %f, %f)\n", b.x, b.y, b.z);
	vector<point> p1, p2; // p1是直线左边所有点集合，p2是直线右边所有点集合
	Point mid = {(a.x+b.x)/2, (a.y+b.y)/2, (a.z+b.z)/2}; // 线段中点
	Point n2 = {b.x-a.x, b.y-a.y, b.z-a.z}; //两个极值点的线段所在的向量
	Point n3 = {n1.y*n2.z-n2.y*n1.z, n2.x*n1.z-n1.x*n2.z, n1.x*n2.y-n2.x*n1.y}; // 计算所在平面内的线段的法向量
	for (int i = 0; i < allPoints.size(); ++i)
	{
		point temp = {allPoints[i].x-mid.x, allPoints[i].y-mid.y, allPoints[i].z-mid.z}; //点集合中任意一个点到直线中点的向量
		double value = n3.x*temp.x + n3.y*temp.y + n3.z*temp.z; //向量和平面内直线法向量的点积
		if(value > threHold) p1.push_back(allPoints[i]);
		else if(value < -threHold) p2.push_back(allPoints[i]);
	}
	printf("left-part:%d, right-part:%d\n", p1.size(), p2.size());
	FindPoint2(p1, a, b, mid, polygon);
	printf("左边的算完了~~~~~");
	FindPoint2(p2, a, b, mid, polygon);
	printf("两边都算完了~~我猜是看不到这一句的把。~~~");
}


/*------《测试篇》-------------我--------是---------分--------割---------线--------------------- */


void drawTest() {
	for (int i=0; i<vtemp.size(); i++) {
		for(int j=i+1; j<vtemp.size(); j++) {
			drawLines(vtemp[i], vtemp[j]);
		}
	}
}

void testSurfs(int i) {
	//for (int i=0; i<vs.size(); i++) {
		printf("%f %f %f %f %f %f\n", vs[i].xmin, vs[i].ymin,  vs[i].zmin, vs[i].xmax, vs[i].ymax,  vs[i].zmax);
		printf(vs[i].isGate ? "true ":"false ");
		printf("%f %f %f %f %f %f\n", vs[i].gxmin, vs[i].gymin,  vs[i].gzmin, vs[i].gxmax, vs[i].gymax,  vs[i].gzmax);
		printf("%f %f %f %f\n", vs[i].a,vs[i].b,vs[i].c,vs[i].d);
		printf("\n");
	//}
}

void testPoints() {
	struct point hp, next;
	for (int i=0; i<vls.size(); i++) {
		printf("\n(%.2f,%.2f,%.2f)-(%.2f,%.2f,%.2f)\n",vls[i].inner[0].x,vls[i].inner[0].y,vls[i].inner[0].z,vls[i].inner[1].x,vls[i].inner[1].y,vls[i].inner[1].z);
		if (!collision(vls[i].inner[0], vls[i].inner[1], &hp, &next, i)) {
			printf("   collision: (%.2f,%.2f,%.2f)\n", hp.x, hp.y, hp.z);
			vtemp.push_back(hp);
		}
	}

}

void testSurfPolygons(vector<point> &allPoints, vector<point> &polygon, point n1) {
	/*for(int i=0; i<allPoints.size(); i++) {
		printf("point%d:(%f, %f, %f)\n", i, allPoints[i].x, allPoints[i].y, allPoints[i].z);
	}*/
}

void testSurfPolygons2() {
	int n;
	vector<point> v, polygon;
	scanf("%d", &n);
	for(int i=0; i<n; i++) {
		point p = {};
		int x, y, z;
		scanf("%d %d %d", &x, &y, &z);
		p.x = x * 1.0;
		p.y = y * 1.0;
		p.z = z * 1.0;
		v.push_back(p);
	}
	point n1 = {0, 0, 1.0};
	dealWith(v, vtemp, n1);
	for(int i=0; i<polygon.size(); i++) {
		//printf("point%d:(%f, %f, %f)\n", i, polygon[i].x, polygon[i].y, polygon[i].z);
		printf("point%d:(%f, %f, %f)\n", i, vtemp[i].x, vtemp[i].y, vtemp[i].z);
	}
}
