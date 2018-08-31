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
vector<surf> vs; // ��������ƽ�������
vector<lines> vls;
double moveLen = 1.0, width=2.0, swidth=0.5;
bool decrease = true;
double percent = 0.9;


// ������
int main(int argc,char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(1200,800);
	glutInitWindowPosition(100,0);
	glutCreateWindow("RFID-model views");
	init();

	dataPrepare(); //׼��һЩ����

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

//����ʾ����
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

	//��ʾ�����ἰ��Χ����
	drawRefections();	
	drawEnvironments();	
	//drawRefections();	
	//drawTest();
	glutSwapBuffers();
}

//����һ��ʼ��������������������뾶��Ϣ
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

//��ʼ������
void init(void)
{
	glEnable(GL_MAP2_VERTEX_3);
	glEnable(GL_DEPTH_TEST);
	//glColor4f(1.0f, 1.0f, 1.0f, 0.5f);//��ɫ0.5 alphaֵ
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);//��Ϻ���
	glEnable(GL_BLEND);
}

//��׽����ƶ�����
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

//�ػ��ƺ���
void reshape(int w, int h)
{
	glViewport(0, 0, (GLint)w, (GLint)h);
}

//�������ú���
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

//����ƶ�����
void mouse(int button, int state, int x, int y)
{
	mx = (GLint)x;
	my = (GLint)y;
}

//���������ắ��
void drawCoordinates(void)
{
	glColor3f(1.0f,0.0f,0.0f); //����ɫ��x��
	glBegin(GL_LINES);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(1.0f, 0.0f, 0.0f);
	glEnd();
	glColor3f(0.0,1.0,0.0); //����ɫ��y��
	glBegin(GL_LINES);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 1.0f, 0.0f);
	glEnd();
	glColor3f(0.0,0.0,1.0); //����ɫ��z��
	glBegin(GL_LINES);
		glVertex3f(0.0f, 0.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, 1.0f);
	glEnd();

}

//��������ֱ�ߺ���
void drawLines(point p1, point p2) {
	glColor3f( 1.0, 1.0, 0.0 );
	glLineWidth(lineWidth);
	glBegin(GL_LINES);
		glVertex3f( p1.x, p1.y, p1.z );
		glVertex3f(  p2.x, p2.y, p2.z );
	glEnd();
}

//����һ���㺯��
void drawPoint(point p) {
	glPointSize(5.0);
	glColor3f( 1.0, 1.0, 0.0 );
	glBegin(GL_POINTS);
		glVertex3f(p.x, p.y, p.z);
	glEnd();
}

//��ʼ������ƽ����Ϣ��
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

// ���ƻ���ƽ����Ϣ
void drawSurfaces() {
	
	glColor4f( 0.5f, 0.5f, 0.5f,0.5f );
	//���Ʊ�����
	glBegin( GL_QUADS );
		glVertex3f( width+moveLen, width/2, -width/2 );
		glVertex3f( width+moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, -width/2 );
	glEnd();
	glColor4f( 0.8f, 0.8f, 0.8f,0.5f );
	//��������
	glBegin( GL_QUADS );
		glVertex3f( moveLen, width/2, -width/2 );
		glVertex3f( moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, width/2, -width/2 );
	glEnd();
	//��������
	glBegin( GL_QUADS );
		glVertex3f( moveLen, -width/2, -width/2 );
		glVertex3f( moveLen, -width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, -width/2 );
	glEnd();
	////����ǰ��
	glBegin( GL_QUADS );
		glVertex3f( moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, width/2, width/2 );
		glVertex3f( width+moveLen, -width/2, width/2 );
		glVertex3f( moveLen, -width/2, width/2 );
	glEnd();
	////���Ʊ�����
	glBegin( GL_QUADS );
		glVertex3f( moveLen, width/2, -width/2 );
		glVertex3f( width+moveLen, width/2, -width/2 );
		glVertex3f( width+moveLen, -width/2, -width/2 );
		glVertex3f( moveLen, -width/2, -width/2 );
	glEnd();

	//��������
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

//����������߿�ֱ��
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


		//�ڿ�
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

//���ƻ�����Ϣ
void drawEnvironments() {
	drawCoordinates(); // ����������
	drawSurfaces(); // ����ƽ����Ϣ
	drawLines(); // ����ƽ�������ֱ��
}

//���ڻ��������ߵĺ���
void drawRefections() {
	glColor3f( 1.0, 1.0, 0.0 );
	glLineWidth(0.01);
	for (int i=0; i<vls.size(); i++) {
		for(int j=0; j<vls[i].inner.size()-1; j++) {
			drawLines(vls[i].inner[j], vls[i].inner[j+1]);
		}

	}
}

// �����������������ߣ�
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

//����׼������
void dataPrepare(){
	//createSphereArray(20, 6); //������������е�����
	initSurfs(); //�����ʼ����ƽ����Ϣ
	//initPoints(); //��ѡ����ֱ�ߵĵ���Ϣ
	initMultipleLight();
	//initSingleLight();
	//testSurfs();
	//testPoints();
	analyze();
}

//���������ߵĺ���
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

//�������������ߵķ���
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
* ��֤�����Ƿ���ĳƽ��֮�ڵĺ�����
**/
int validate(point cross, surf temp) {
	//printf("�ͻ�����ײ���������ǣ�(%f,%f,%f)\n", cross.x,cross.y,cross.z);
	//printf("��ײ����Ϣx��[%.2f %.2f]��y��[%.2f %.2f]��z��[%.2f %.2f]\n", temp.xmin,temp.xmax,temp.ymin,temp.ymax,temp.zmin,temp.zmax);
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

//����p1��p2������s֮��Ľ������꺯�������ؽ��� temp
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
*  ����ײ����0�����ѽ������Ϣ������hp���棬�������Ϣ������next����򷵻�1��
*  p1�ǵ�һ���㣬p2�ǵڶ����㣬hp�ǽ���㣬next����һ������㣬
*/
int collision(point p1, point p2, struct point *hp, struct point *next, int index) {
	for (int i=0; i<vs.size(); i++) {
		struct surf temp = vs[i];
		double result = (temp.a * p1.x + temp.b * p1.y + temp.c * p1.z + temp.d) * (temp.a * p2.x + temp.b * p2.y + temp.c * p2.z + temp.d);
		
		if (result <= 0) {
			//printf("��ײ�ˡ�����\n");
			int size = vls[index].surfIndex.size();
			if (result == 0) {
				if (size >=1 && vls[index].surfIndex[size-1] == i) continue;
				else if (size >=2 && vls[index].surfIndex[size-2] == i) continue;
			} 
			struct point cross = getCross(p1, p2, temp);

			//printf("�ͻ�����ײ���������ǣ�(%.2f,%.2f,%.2f)\n", cross.x,cross.y,cross.z);
			//printf("�ͻ�����ײ���������ǣ�(%f,%f,%f)\n", cross.x,cross.y,cross.z);
			
			//printf("��ײ����Ϣx��[%.2f %.2f]��y��[%.2f %.2f]��z��[%.2f %.2f]\n", temp.xmin,temp.xmax,temp.ymin,temp.ymax,temp.zmin,temp.zmax);
			if (validate(cross, temp)) {
					//printf("��ײ����������ײ����ܷ�Χ����,������ײ����");
					//printf(temp.isGate ? "��������\n":"��ͨ����\n");
					if (temp.isGate && cross.x >= temp.gxmin && cross.x <= temp.gxmax && 
						cross.y >= temp.gymin && cross.y <= temp.gymax &&
						cross.z >= temp.gzmin && cross.z <= temp.gzmax) continue;
				//printf("��ײ����������ײ�����Ч��Χ����\n");
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
				//printf("������ǣ�%f %f %f��,��ײ���ǵ�%d��ƽ��\n", next->x, next->y,next->z,i);
				vls[index].surfIndex.push_back(i);
				return 0;
			}			
		}
	}
	return 1;
}

//ѭ�����亯��
void analyze() {
	int c = 3;
	for (int i=0; i<vls.size(); i++) {
		//printf("��%d��ֱ��...\n", i);
		int index = 0;
		struct point hp, next;
		struct lines line = vls[i];
		while (index < vls[i].inner.size()-1) {
			//if (index > c) break;
			//printf("�ڼ����㣺%d\n", index);
			/*printf("ֱ����ʼ����ֹ����ֱ��ǣ�%.2f, %.2f, %.2f������%.2f, %.2f, %.2f��\n", vls[i].inner[index].x, vls[i].inner[index].y, vls[i].inner[index].z,
					vls[i].inner[index+1].x, vls[i].inner[index+1].y, vls[i].inner[index+1].z);*/
			if (!collision(vls[i].inner[index], vls[i].inner[index+1], &hp, &next, i)) {
				//printf("��%d��������ǣ�%f, %f, %f��\n", index, hp.x, hp.y, hp.z);
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


/*------������ƪ��-------------��--------��---------��--------��---------��--------------------- */


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

