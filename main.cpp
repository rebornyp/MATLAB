#include <stdio.h>
#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <glut.h>
#include <vector>
using namespace std;

#define pi 3.1415926
#define threHold 0.001

struct point {
	double x, y, z;
};

struct surf {
	double xmin, xmax, ymin, ymax, zmin, zmax;
	double a,b,c,d;
	bool isGate;
	double gxmin, gxmax, gymin, gymax, gzmin, gzmax;
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
void testSurfs();
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

GLint mx,my; //position of mouse
GLint m_state=0; //mouse usage
GLfloat x_angle=0.0f, y_angle=0.0f; //angle of eye
GLfloat dist=10.0f; //distance from the eye
double lineWidth = 0.01;
vector<point> vp, vtemp;
vector<surf> vs; // 保存所有平面的数组
vector<lines> vls;
double moveLen = 1.0, width=2.0, swidth=0.5;
bool decrease = true;
double percent = 0.9;


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
	glColor3f( 1.0, 1.0, 0.0 );
	glLineWidth(lineWidth);
	glBegin(GL_LINES);
		glVertex3f( p1.x, p1.y, p1.z );
		glVertex3f(  p2.x, p2.y, p2.z );
	glEnd();
}

//绘制一个点函数
void drawPoint(point p) {
	glPointSize(5.0);
	glColor3f( 1.0, 1.0, 0.0 );
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
}

//后期绘制折射线的函数
void drawRefections() {
	glColor3f( 1.0, 1.0, 0.0 );
	glLineWidth(0.01);
	for (int i=0; i<vls.size(); i++) {
		for(int j=0; j<vls[i].inner.size()-1; j++) {
			drawLines(vls[i].inner[j], vls[i].inner[j+1]);
		}

	}
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
	//initPoints(); //节选入射直线的点信息
	initMultipleLight();
	//initSingleLight();
	//testSurfs();
	//testPoints();
	analyze();
}

//测试入射线的函数
void initMultipleLight() {
	double k = 6;
	struct point p1={0,0.3,0};
	int n = 10;
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
	double k = 10;
	struct point p1={0,0.3,0};
	struct point p2={(k+1)*moveLen-k*p1.x, (k+1)*-0.15-k*p1.y, (k+1)*0-k*p1.z};
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
*  p1是第一个点，p2是第二个点，hp是交叉点，next是下一个折射点，
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


/*------《测试篇》-------------我--------是---------分--------割---------线--------------------- */


void drawTest() {
	for (int i=0; i<vtemp.size(); i++) {
		drawPoint(vtemp[i]);
	}
}

void testSurfs() {
	for (int i=0; i<vs.size(); i++) {
		printf("%f %f %f %f %f %f\n", vs[i].xmin, vs[i].ymin,  vs[i].zmin, vs[i].xmax, vs[i].ymax,  vs[i].zmax);
		printf(vs[i].isGate ? "true ":"false ");
		printf("%f %f %f %f %f %f\n", vs[i].gxmin, vs[i].gymin,  vs[i].gzmin, vs[i].gxmax, vs[i].gymax,  vs[i].gzmax);
		printf("%f %f %f %f\n", vs[i].a,vs[i].b,vs[i].c,vs[i].d);
		printf("\n");
	}
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

