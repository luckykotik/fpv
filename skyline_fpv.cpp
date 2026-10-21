// SKYLINE FPV — симулятор FPV-дрона. C++17, SDL2 + OpenGL (фиксированный конвейер).
// Linux:   g++ -O2 -std=c++17 skyline_fpv.cpp -o skyline $(sdl2-config --cflags --libs) -lGL
// macOS:   clang++ -O2 -std=c++17 skyline_fpv.cpp -o skyline $(sdl2-config --cflags --libs) -framework OpenGL
// Windows (MSYS2 mingw64): g++ -O2 -std=c++17 skyline_fpv.cpp -o skyline.exe $(sdl2-config --cflags --libs) -lopengl32
#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cctype>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;
typedef double D;
static const D PI=3.14159265358979323846, G=9.81;

// ---------- математика ----------
struct V{D x,y,z;};
V operator+(V a,V b){return{a.x+b.x,a.y+b.y,a.z+b.z};}
V operator-(V a,V b){return{a.x-b.x,a.y-b.y,a.z-b.z};}
V operator*(V a,D s){return{a.x*s,a.y*s,a.z*s};}
D dot(V a,V b){return a.x*b.x+a.y*b.y+a.z*b.z;}
V cross(V a,V b){return{a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
D len(V a){return sqrt(dot(a,a));}
V norm(V a){D l=len(a);return l>1e-9?a*(1/l):V{0,1,0};}
struct Q{D w,x,y,z;};
Q mul(Q a,Q b){return{a.w*b.w-a.x*b.x-a.y*b.y-a.z*b.z,a.w*b.x+a.x*b.w+a.y*b.z-a.z*b.y,a.w*b.y-a.x*b.z+a.y*b.w+a.z*b.x,a.w*b.z+a.x*b.y-a.y*b.x+a.z*b.w};}
V rot(Q q,V v){V u{q.x,q.y,q.z};V t=cross(u,v)*2;return v+t*q.w+cross(u,t);}
Q axisQ(V a,D ang){D s=sin(ang/2);return{cos(ang/2),a.x*s,a.y*s,a.z*s};}
Q qnorm(Q q){D l=sqrt(q.w*q.w+q.x*q.x+q.y*q.y+q.z*q.z);return{q.w/l,q.x/l,q.y/l,q.z/l};}
Q qconj(Q q){return{q.w,-q.x,-q.y,-q.z};}
static D R01(){return rand()/(D)RAND_MAX;}
static D clampd(D x,D a,D b){return x<a?a:(x>b?b:x);}
string fmt(const char*f,...){char b[96];va_list a;va_start(a,f);vsnprintf(b,96,f,a);va_end(a);return b;}

// ---------- мир ----------
D hgt(D x,D z){D d=hypot(x,z),h=sin(x*.02)*6+cos(z*.017)*7+sin((x+z)*.045)*3+sin(x*.09)*cos(z*.08)*1.2;
  return h*clampd((d-20)/50,0,1)+sin(x*.006+1)*cos(z*.007)*14*clampd((d-60)/100,0,1);}
struct Batch{vector<float> v,c;int n=0;
  void vtx(V p,float r,float g,float b){v.push_back(p.x);v.push_back(p.y);v.push_back(p.z);c.push_back(r);c.push_back(g);c.push_back(b);n++;}
  void tri(V a,V b,V d,float r,float g,float bl){vtx(a,r,g,bl);vtx(b,r,g,bl);vtx(d,r,g,bl);}
  void quad(V a,V b,V c2,V d,float r,float g,float bl){tri(a,b,c2,r,g,bl);tri(a,c2,d,r,g,bl);}
  void draw(){glEnableClientState(GL_VERTEX_ARRAY);glEnableClientState(GL_COLOR_ARRAY);glVertexPointer(3,GL_FLOAT,0,v.data());glColorPointer(3,GL_FLOAT,0,c.data());glDrawArrays(GL_TRIANGLES,0,n);glDisableClientState(GL_COLOR_ARRAY);glDisableClientState(GL_VERTEX_ARRAY);}
};
Batch terr,props,sky;
struct Bld{D x,z,w,d,top;};struct Tr{D x,z,r,top;};struct Gate{D x,y,z,yaw;};
vector<Bld> blds;vector<Tr> trees;vector<Gate> gates;
static const V Ld={.4,.75,.35};

void buildTerrain(){
  const int N=160;const D C=5,O=-400;static D hh[N+1][N+1];V L=norm(Ld);
  for(int i=0;i<=N;i++)for(int j=0;j<=N;j++)hh[i][j]=hgt(O+i*C,O+j*C);
  auto vert=[&](int i,int j){D x=O+i*C,z=O+j*C,h=hh[i][j];
    D hx=hh[min(i+1,N)][j]-hh[max(i-1,0)][j],hz=hh[i][min(j+1,N)]-hh[i][max(j-1,0)];
    V n=norm({-hx,2*C,-hz});float l=.5f+.5f*(float)max(0.,dot(n,L));
    D nz=sin(x*.31)*cos(z*.27)*.5+.5;float r,g,b;
    if(h<-2.2){r=.84f;g=.77f;b=.56f;}else if(h<9){r=.2f+.07f*nz;g=.42f+.13f*nz;b=.14f;}else if(h<16){r=.5f;g=.46f+.04f*nz;b=.4f;}else{r=.93f;g=.95f;b=.97f;}
    terr.vtx({x,h,z},r*l,g*l,b*l);};
  for(int i=0;i<N;i++)for(int j=0;j<N;j++){vert(i,j);vert(i,j+1);vert(i+1,j);vert(i+1,j);vert(i,j+1);vert(i+1,j+1);}
}
void cone(D x,D y,D z,D r,D h,float R,float Gc,float B){const int s=7;
  for(int k=0;k<s;k++){D a=k*2*PI/s,b=(k+1)*2*PI/s;float sh=.7f+.3f*(float)cos(a-1);
    props.tri({x,y+h,z},{x+r*cos(a),y,z+r*sin(a)},{x+r*cos(b),y,z+r*sin(b)},R*sh,Gc*sh,B*sh);}}
void tree(D x,D y,D z,D s){
  props.quad({x-.2*s,y,z},{x+.2*s,y,z},{x+.2*s,y+2*s,z},{x-.2*s,y+2*s,z},.35f,.24f,.14f);
  props.quad({x,y,z-.2*s},{x,y,z+.2*s},{x,y+2*s,z+.2*s},{x,y+2*s,z-.2*s},.3f,.2f,.12f);
  cone(x,y+1.5*s,z,2.2*s,4*s,.1f,.33f,.15f);cone(x,y+3.5*s,z,1.6*s,3.5*s,.12f,.38f,.17f);}
struct Face{V o,t;D L;V n;};
void building(D cx,D cz,D w,D d,D y0,D y1){
  Face fs[4]={{{cx-w/2,y0,cz+d/2},{1,0,0},w,{0,0,1}},{{cx+w/2,y0,cz+d/2},{0,0,-1},d,{1,0,0}},
              {{cx+w/2,y0,cz-d/2},{-1,0,0},w,{0,0,-1}},{{cx-w/2,y0,cz-d/2},{0,0,1},d,{-1,0,0}}};
  V up={0,y1-y0,0};
  for(int i=0;i<4;i++){Face&f=fs[i];float sh=.6f+.12f*i;
    props.quad(f.o,f.o+f.t*f.L,f.o+f.t*f.L+up,f.o+up,.55f*sh,.6f*sh,.68f*sh);
    for(D y=y0+3;y<y1-2;y+=3)for(D s=1.2;s<f.L-2.4;s+=2.4){bool on=R01()<.4;V o=f.o+f.t*s+f.n*.05+V{0,y,0};
      props.quad(o,o+f.t*1.2,o+f.t*1.2+V{0,1.5,0},o+V{0,1.5,0},on?1.f:.1f,on?.88f:.16f,on?.55f:.26f);}}
  props.quad({cx-w/2,y1,cz-d/2},{cx+w/2,y1,cz-d/2},{cx+w/2,y1,cz+d/2},{cx-w/2,y1,cz+d/2},.35f,.35f,.38f);
  blds.push_back({cx,cz,w,d,y1});}
void buildProps(){
  for(D bx:{100.,125.,150.})for(D bz:{-110.,-85.,-60.}){D w=9+R01()*3,d=9+R01()*3,hh=16+R01()*24;
    D base=min(min(hgt(bx-w/2,bz-d/2),hgt(bx+w/2,bz+d/2)),min(hgt(bx-w/2,bz+d/2),hgt(bx+w/2,bz-d/2)))-2;building(bx,bz,w,d,base,base+hh);}
  for(int n=0,k=0;n<520&&k<9000;k++){D x=(R01()-.5)*780,z=(R01()-.5)*780,h=hgt(x,z);
    if(h<-1||h>13||hypot(x,z)<24||(x>80&&x<180&&z>-130&&z<-30))continue;
    D s=.8+R01()*1.1;tree(x,h-.3,z,s);trees.push_back({x,z,.8*s,h+7*s});n++;}
  D gw[10][3]={{0,-30,5},{25,-70,8},{65,-85,12},{112,-72,10},{137,-97,12},{175,-75,14},{170,0,16},{110,60,12},{40,70,10},{-20,40,8}};
  for(int i=0;i<10;i++){D*a=gw[i],*b=gw[(i+1)%10];gates.push_back({a[0],hgt(a[0],a[1])+a[2],a[1],atan2(b[0]-a[0],b[1]-a[1])});}
  // небо
  const int RN=10,S=28;
  for(int i=-1;i<RN;i++){D a0=i<0?-.12:i*PI/2/RN,a1=(i+1)*PI/2/RN;
    for(int k=0;k<S;k++){D p0=k*2*PI/S,p1=(k+1)*2*PI/S;
      auto P=[&](D a,D b){return V{900*cos(a)*cos(b),900*sin(a),900*cos(a)*sin(b)};};
      auto C=[&](D a,float*c){D t=pow(max(0.,sin(a)),.5);c[0]=(float)(.78*(1-t)+.12*t);c[1]=(float)(.87*(1-t)+.35*t);c[2]=(float)(.95*(1-t)+.78*t);};
      float c0[3],c1[3];C(a0,c0);C(a1,c1);
      sky.vtx(P(a0,p0),c0[0],c0[1],c0[2]);sky.vtx(P(a0,p1),c0[0],c0[1],c0[2]);sky.vtx(P(a1,p0),c1[0],c1[1],c1[2]);
      sky.vtx(P(a0,p1),c0[0],c0[1],c0[2]);sky.vtx(P(a1,p1),c1[0],c1[1],c1[2]);sky.vtx(P(a1,p0),c1[0],c1[1],c1[2]);}}
}

// ---------- шрифт 5x7 (ASCII 32..90) ----------
static const unsigned char FONT[59][5]={
{0,0,0,0,0},{0,0,0x5F,0,0},{0,7,0,7,0},{0x14,0x7F,0x14,0x7F,0x14},{0x24,0x2A,0x7F,0x2A,0x12},{0x23,0x13,8,0x64,0x62},{0x36,0x49,0x55,0x22,0x50},{0,5,3,0,0},
{0,0x1C,0x22,0x41,0},{0,0x41,0x22,0x1C,0},{0x14,8,0x3E,8,0x14},{8,8,0x3E,8,8},{0,0x50,0x30,0,0},{8,8,8,8,8},{0,0x60,0x60,0,0},{0x20,0x10,8,4,2},
{0x3E,0x51,0x49,0x45,0x3E},{0,0x42,0x7F,0x40,0},{0x42,0x61,0x51,0x49,0x46},{0x21,0x41,0x45,0x4B,0x31},{0x18,0x14,0x12,0x7F,0x10},{0x27,0x45,0x45,0x45,0x39},{0x3C,0x4A,0x49,0x49,0x30},{1,0x71,9,5,3},{0x36,0x49,0x49,0x49,0x36},{6,0x49,0x49,0x29,0x1E},
{0,0x36,0x36,0,0},{0,0x56,0x36,0,0},{8,0x14,0x22,0x41,0},{0x14,0x14,0x14,0x14,0x14},{0,0x41,0x22,0x14,8},{2,1,0x51,9,6},{0x32,0x49,0x79,0x41,0x3E},
{0x7E,0x11,0x11,0x11,0x7E},{0x7F,0x49,0x49,0x49,0x36},{0x3E,0x41,0x41,0x41,0x22},{0x7F,0x41,0x41,0x22,0x1C},{0x7F,0x49,0x49,0x49,0x41},{0x7F,9,9,9,1},{0x3E,0x41,0x49,0x49,0x7A},
{0x7F,8,8,8,0x7F},{0,0x41,0x7F,0x41,0},{0x20,0x40,0x41,0x3F,1},{0x7F,8,0x14,0x22,0x41},{0x7F,0x40,0x40,0x40,0x40},{0x7F,2,0xC,2,0x7F},{0x7F,4,8,0x10,0x7F},{0x3E,0x41,0x41,0x41,0x3E},
{0x7F,9,9,9,6},{0x3E,0x41,0x51,0x21,0x5E},{0x7F,9,0x19,0x29,0x46},{0x46,0x49,0x49,0x49,0x31},{1,1,0x7F,1,1},{0x3F,0x40,0x40,0x40,0x3F},{0x1F,0x20,0x40,0x20,0x1F},{0x3F,0x40,0x38,0x40,0x3F},
{0x63,0x14,8,0x14,0x63},{7,8,0x70,8,7},{0x61,0x51,0x49,0x45,0x43}};
void rect2(D x0,D y0,D x1,D y1){glVertex2d(x0,y0);glVertex2d(x1,y0);glVertex2d(x1,y1);glVertex2d(x0,y1);}
void txt(D x,D y,const string&s,D sc,float r,float g,float b,float a=1){
  glColor4f(r,g,b,a);glBegin(GL_QUADS);
  for(char ch:s){int c=toupper((unsigned char)ch)-32;
    if(c>=0&&c<59)for(int i=0;i<5;i++)for(int j=0;j<7;j++)if((FONT[c][i]>>j)&1)rect2(x+i*sc,y+j*sc,x+(i+1)*sc,y+(j+1)*sc);
    x+=6*sc;}
  glEnd();}
void txs(D x,D y,const string&s,D sc,float r=.78f,float g=1,float b=.78f){txt(x+sc*.8,y+sc*.8,s,sc,0,0,0,.75f);txt(x,y,s,sc,r,g,b);}
D tw(const string&s,D sc){return s.size()*6*sc;}

// ---------- состояние дрона ----------
V p={0,.06,0},v={0,0,0};Q q={1,0,0,0};
D thrust=0,tm=0,mah=0,curI=0,flightT=0,crashT=0,msgT=0,kthr=0,aRpm=0,prevDd=0,finishT=0;
bool armed=false,onGround=true,paused=false,launched=false,finished=false,gpOn=false;
int mode=0,camM=0,gi=0;string msg;
static const D MAXA=2.2*G,CAPMAH=1500,TILT=30*PI/180;
struct In{D thr=0,yaw=0,pit=0,rol=0;}in;
SDL_Window*win;SDL_GameController*pad=nullptr;

void say(const string&s,D t=2){msg=s;msgT=t;}
void resetAll(bool full){p={0,.06,0};v={0,0,0};q={1,0,0,0};thrust=0;armed=false;onGround=true;crashT=0;kthr=mode?-1:0;
  if(full){mah=0;flightT=0;gi=0;launched=false;finished=false;prevDd=0;}}
void crash(const char*w){if(crashT>0)return;crashT=2.4;armed=false;v={0,0,0};thrust=0;say(string("CRASH: ")+w,2.2);}
D voltage(){D rem=max(0.,1-mah/CAPMAH);return 13.2+3.6*pow(rem,.6)-curI*.012;}

D ex(D x){return x*(.35+.65*fabs(x));}
D ramp(D a,D b,D dt){D r=b==0?9:5;return a+clampd(b-a,-r*dt,r*dt);}
D axn(SDL_GameControllerAxis a){D x=SDL_GameControllerGetAxis(pad,a)/32767.;if(fabs(x)<.08)return 0;return clampd((x-(x>0?.08:-.08))/.92,-1,1);}
void gather(D dt){
  const Uint8*k=SDL_GetKeyboardState(0);
  auto ax=[&](int a,int b){return (k[a]?1.:0.)-(k[b]?1.:0.);};
  D tt=ax(SDL_SCANCODE_UP,SDL_SCANCODE_DOWN),yy=ax(SDL_SCANCODE_RIGHT,SDL_SCANCODE_LEFT)+ax(SDL_SCANCODE_E,SDL_SCANCODE_Q);
  D pp=ax(SDL_SCANCODE_W,SDL_SCANCODE_S),rr=ax(SDL_SCANCODE_D,SDL_SCANCODE_A);
  if(tt!=0)gpOn=false;
  if(mode==0)kthr=ramp(kthr,tt,dt);else kthr=clampd(kthr+tt*.8*dt,-1,1);
  static D y_=0,p_=0,r_=0;y_=ramp(y_,yy,dt);p_=ramp(p_,pp,dt);r_=ramp(r_,rr,dt);
  in.thr=kthr;in.yaw=y_;in.pit=p_;in.rol=r_;
  if(pad){D lx=axn(SDL_CONTROLLER_AXIS_LEFTX),ly=-axn(SDL_CONTROLLER_AXIS_LEFTY),rx=axn(SDL_CONTROLLER_AXIS_RIGHTX),ry=-axn(SDL_CONTROLLER_AXIS_RIGHTY);
    if(fabs(ly)>.2||fabs(lx)>.2||fabs(rx)>.2||fabs(ry)>.2)gpOn=true;
    if(gpOn)in.thr=ly;
    in.yaw=clampd(in.yaw+lx,-1,1);in.rol=clampd(in.rol+rx,-1,1);in.pit=clampd(in.pit+ry,-1,1);}
}
void toggleArm(){
  if(armed){armed=false;say("DISARMED");return;}
  if(crashT>0)return;
  if(mode==1&&in.thr>-.6){say("LOWER THROTTLE TO ARM",2);return;}
  armed=true;say("ARMED",1.5);}
void setMode(){mode^=1;kthr=mode?-1:0;say(mode?"ACRO MODE":"ANGLE MODE",1.5);}

// ---------- физика ----------
void step(D dt){
  tm+=dt;if(msgT>0)msgT-=dt;
  if(crashT>0){crashT-=dt;if(crashT<=0)resetAll(false);return;}
  V wind={sin(tm*.11)*2.5+sin(tm*.9)*.8,0,cos(tm*.07)*2+sin(tm*.6)*.7};
  V up=rot(q,{0,1,0});D rem=1-mah/CAPMAH;
  D cmd=0;
  if(armed){
    if(mode==0){D vt=in.thr*(in.thr>0?5:3);D need=(G+3*(vt-v.y))/max(.45,up.y);cmd=clampd(need/MAXA,0,1);}
    else cmd=clampd((in.thr+1)/2,0,1);
    if(rem<=0)cmd=min(cmd,.3);}
  thrust+=(cmd-thrust)*(1-exp(-dt/.05));
  V w={0,0,0};
  if(onGround&&thrust<.3){V f=rot(q,{0,0,-1});q=axisQ({0,1,0},atan2(-f.x,-f.z));}
  else{
    if(mode==1){w={-ex(in.pit)*10.5,-ex(in.yaw)*7,-ex(in.rol)*10.5};}
    else{V f=rot(q,{0,0,-1});V fh=norm({f.x,0,f.z});if(hypot(f.x,f.z)<.05)fh=norm(rot(q,{0,1,0}));
      V rh={-fh.z,0,fh.x};V d=norm(V{0,1,0}+fh*tan(ex(in.pit)*.75)+rh*tan(ex(in.rol)*.75));
      V eb=rot(qconj(q),cross(up,d));w={clampd(eb.x*9,-10,10),-ex(in.yaw)*3.5,clampd(eb.z*9,-10,10)};}
    w=w+V{R01()-.5,R01()-.5,R01()-.5}*(.12*thrust);}
  Q dq=mul(q,Q{0,w.x,w.y,w.z});q=qnorm(Q{q.w+dq.w*.5*dt,q.x+dq.x*.5*dt,q.y+dq.y*.5*dt,q.z+dq.z*.5*dt});
  up=rot(q,{0,1,0});
  V rel=v-wind;V acc=up*(thrust*MAXA)+V{0,-G,0}-rel*(.25+.02*len(rel));
  v=v+acc*dt;p=p+v*dt;
  curI=armed?1+45*pow(thrust,1.5):0;mah+=curI*dt/3.6;if(armed&&!onGround)launched=true;
  if(armed&&launched&&!finished)flightT+=dt;
  // столкновения
  D g=hgt(p.x,p.z);
  if(p.y<max(g,-3.)+.06){
    if(g<-3)crash("WATER");
    else if(!onGround&&(-v.y>6.5||hypot(v.x,v.z)>12||up.y<.35))crash("GROUND");
    else{p.y=g+.06;if(v.y<0)v.y=0;onGround=true;D f=exp(-dt*5);v.x*=f;v.z*=f;}
  }else onGround=false;
  for(auto&b:blds)if(fabs(p.x-b.x)<b.w/2+.3&&fabs(p.z-b.z)<b.d/2+.3&&p.y<b.top+.2){crash("BUILDING");break;}
  for(auto&t:trees){D dx=p.x-t.x,dz=p.z-t.z;if(dx*dx+dz*dz<t.r*t.r+.2&&p.y<t.top){crash("TREE");break;}}
  if(fabs(p.x)>395||fabs(p.z)>395){p.x=clampd(p.x,-395,395);p.z=clampd(p.z,-395,395);say("OUT OF RANGE",1);}
  // ворота
  if(!finished&&gi<(int)gates.size()){Gate&gt=gates[gi];V n={sin(gt.yaw),0,cos(gt.yaw)};V c={gt.x,gt.y,gt.z};
    D dd=dot(p-c,n),lat=dot(p-c,V{n.z,0,-n.x}),vt=p.y-gt.y;
    if(prevDd*dd<0&&fabs(lat)<2.9&&fabs(vt)<2.9){gi++;if(gi>=(int)gates.size()){finished=true;finishT=flightT;say(fmt("FINISH %d:%02d",(int)finishT/60,(int)finishT%60),6);}else say(fmt("GATE %d/%d",gi,(int)gates.size()),1);
      if(gi<(int)gates.size()){Gate&g2=gates[gi];prevDd=dot(p-V{g2.x,g2.y,g2.z},V{sin(g2.yaw),0,cos(g2.yaw)});}}
    else prevDd=dd;}
}

// ---------- рендер мира ----------
void setView(V r,V u,V f,V pos){D m[16]={r.x,u.x,-f.x,0,r.y,u.y,-f.y,0,r.z,u.z,-f.z,0,0,0,0,1};glMultMatrixd(m);glTranslated(-pos.x,-pos.y,-pos.z);}
void persp(D fov,D asp){D t=.1*tan(fov*PI/360),r=t*asp;glFrustum(-r,r,-t,t,.1,1500);}
void disc(D x,D y,D z,D r,float a){glColor4f(0,0,0,a);glBegin(GL_TRIANGLE_FAN);glVertex3d(x,y,z);for(int i=0;i<=16;i++)glVertex3d(x+r*cos(i*PI/8),y,z+r*sin(i*PI/8));glEnd();}
void drawGate(const Gate&g,bool act){
  glPushMatrix();glTranslated(g.x,g.y,g.z);glRotated(g.yaw*180/PI,0,1,0);
  float a=(float)(act?1:.35);if(act)glColor4f(1,.55f+.3f*(float)sin(tm*8),.05f,1);else glColor4f(.2f,.9f,1,a);
  const D h=2.8,w=.35;glBegin(GL_QUADS);
  auto R=[&](D x0,D y0,D x1,D y1){glVertex3d(x0,y0,0);glVertex3d(x1,y0,0);glVertex3d(x1,y1,0);glVertex3d(x0,y1,0);};
  R(-h-w,-h-w,h+w,-h);R(-h-w,h,h+w,h+w);R(-h-w,-h,-h,h);R(h,-h,h+w,h);glEnd();glPopMatrix();}
// ---------- модель дрона: 5" FPV (рама X, моторы, стек, камера, АКБ, пропеллеры) ----------
struct DX{D yaw,pit;V c;};
static V dxf(const DX&t,V l){D cp=cos(t.pit),sp=sin(t.pit),cy=cos(t.yaw),sy=sin(t.yaw);
  D y1=l.y*cp-l.z*sp,z1=l.y*sp+l.z*cp;return{t.c.x+l.x*cy+z1*sy,t.c.y+y1,t.c.z-l.x*sy+z1*cy};}
static V dnf(const DX&t,V l){DX u=t;u.c={0,0,0};return dxf(u,l);}
static float dsh(V n){V w=rot(q,norm(n));return (float)min(1.,.34+.66*max(0.,dot(w,norm(Ld)))+.12*max(0.,w.y));}
void dbox(V c,V s,D yaw,D pit,float R,float Gc,float B,float emis=0){
  DX t{yaw,pit,c};D h[3]={s.x/2,s.y/2,s.z/2};static const int cu[4]={-1,1,1,-1},cw[4]={-1,-1,1,1};
  glBegin(GL_QUADS);
  for(int a=0;a<3;a++)for(int sg=-1;sg<=1;sg+=2){
    D nn[3]={0,0,0};nn[a]=sg;float k=min(1.f,dsh(dnf(t,{nn[0],nn[1],nn[2]}))+emis);glColor3f(R*k,Gc*k,B*k);
    int u=(a+1)%3,w=(a+2)%3;
    for(int i=0;i<4;i++){D pp[3];pp[a]=sg*h[a];pp[u]=cu[i]*h[u];pp[w]=cw[i]*h[w];V o=dxf(t,{pp[0],pp[1],pp[2]});glVertex3d(o.x,o.y,o.z);}}
  glEnd();}
void dcyl(V c,D r,D h,float R,float Gc,float B,int sg=10){
  glBegin(GL_QUADS);
  for(int i=0;i<sg;i++){D a=i*2*PI/sg,b=(i+1)*2*PI/sg,m=(a+b)/2;float k=dsh({cos(m),.15,sin(m)});glColor3f(R*k,Gc*k,B*k);
    glVertex3d(c.x+r*cos(a),c.y,c.z+r*sin(a));glVertex3d(c.x+r*cos(b),c.y,c.z+r*sin(b));
    glVertex3d(c.x+r*cos(b),c.y+h,c.z+r*sin(b));glVertex3d(c.x+r*cos(a),c.y+h,c.z+r*sin(a));
    k=dsh({0,1,0});glColor3f(R*k,Gc*k,B*k);
    glVertex3d(c.x,c.y+h,c.z);glVertex3d(c.x+r*cos(a),c.y+h,c.z+r*sin(a));glVertex3d(c.x+r*cos(b),c.y+h,c.z+r*sin(b));glVertex3d(c.x,c.y+h,c.z);}
  glEnd();}
void dprop(V c,D spin,float tr,float tg,float tb){
  const D Rp=.135,r0=.012;bool blur=aRpm>.05;float al=blur?.55f:1.f;
  auto off0=[&](D t){return -.006+.008*t;};auto off1=[&](D t){return .006+.010*t;};
  glBegin(GL_QUADS);
  for(int k=0;k<3;k++){D a=spin+k*2*PI/3,ca=cos(a),sa=sin(a);
    for(int seg=0;seg<2;seg++){D t0=seg?.72:0,t1=seg?1:.72;
      if(seg)glColor4f(tr,tg,tb,al);else glColor4f(.1f,.11f,.13f,al);
      D ra=r0+t0*(Rp-r0),rb=r0+t1*(Rp-r0);
      auto P=[&](D rr,D o){return V{c.x+ca*rr-sa*o,c.y,c.z+sa*rr+ca*o};};
      V a1=P(ra,off0(t0)),a2=P(ra,off1(t0)),b2=P(rb,off1(t1)),b1=P(rb,off0(t1));
      glVertex3d(a1.x,a1.y,a1.z);glVertex3d(a2.x,a2.y,a2.z);glVertex3d(b2.x,b2.y,b2.z);glVertex3d(b1.x,b1.y,b1.z);}}
  glEnd();
  if(blur){glBegin(GL_TRIANGLE_FAN);glColor4f(.85f,.9f,1,.12f);glVertex3d(c.x,c.y,c.z);
    for(int i=0;i<=20;i++)glVertex3d(c.x+Rp*cos(i*PI/10),c.y,c.z+Rp*sin(i*PI/10));glEnd();
    glLineWidth(2);glBegin(GL_LINE_LOOP);glColor4f(tr,tg,tb,.5f);for(int i=0;i<20;i++)glVertex3d(c.x+Rp*cos(i*PI/10),c.y,c.z+Rp*sin(i*PI/10));glEnd();}}
void drone(D spin){
  glPushMatrix();glTranslated(p.x,p.y,p.z);V r=rot(q,{1,0,0}),u=rot(q,{0,1,0}),b=rot(q,{0,0,1});
  D m[16]={r.x,r.y,r.z,0,u.x,u.y,u.z,0,b.x,b.y,b.z,0,0,0,0,1};glMultMatrixd(m);glScaled(2,2,2);
  bool on=armed&&((int)(tm*4))%2==0;
  // рама: пластины и стойки
  dbox({0,-.014,.005},{.05,.004,.15},0,0,.1f,.1f,.11f);dbox({0,.016,.005},{.05,.004,.15},0,0,.12f,.12f,.13f);
  for(int sx=-1;sx<=1;sx+=2)for(int sz=-1;sz<=1;sz+=2)dbox({sx*.02,.001,.005+sz*.055},{.006,.03,.006},0,0,.85f,.68f,.22f);
  // стек: ESC / FC / VTX
  dbox({0,-.007,.005},{.036,.003,.036},0,0,.05f,.42f,.2f);dbox({0,.001,.005},{.036,.003,.036},0,0,.1f,.16f,.45f);dbox({0,.009,.005},{.036,.003,.036},0,0,.12f,.12f,.14f);
  // руки, моторы, пропеллеры
  for(int i=0;i<4;i++){D a=PI/4+i*PI/2,mx=.25*cos(a),mz=.25*sin(a);bool fr=mz<0;
    dbox({.135*cos(a),-.002,.135*sin(a)},{.27,.006,.024},-a,0,.09f,.09f,.1f);
    dbox({.06*cos(a),.002,.06*sin(a)},{.05,.003,.012},-a,0,.9f,.45f,.1f);
    dcyl({mx,-.017,mz},.008,.014,.06f,.06f,.06f,8);
    dcyl({mx,-.005,mz},.013,.008,.2f,.2f,.22f);dcyl({mx,.003,mz},.0185,.02,.8f,.14f,.1f);dcyl({mx,.023,mz},.006,.008,.85f,.85f,.9f,6);
    float lk=armed?1.f:.25f;if(fr)dbox({mx*.86,-.004,mz*.86},{.008,.004,.008},0,0,.2f*lk,1.f*lk,.3f*lk,.5f);
    else dbox({mx*.86,-.004,mz*.86},{.008,.004,.008},0,0,1.f*lk,.15f*lk,.1f*lk,.5f);}
  // камера (наклон вверх на TILT)
  DX ct{0,TILT,{0,.007,-.07}};dbox(ct.c,{.028,.028,.026},0,TILT,.07f,.07f,.08f);
  for(int s=-1;s<=1;s+=2)dbox(dxf(ct,{s*.016,0,0}),{.004,.03,.03},0,TILT,.9f,.45f,.1f);
  dbox(dxf(ct,{0,0,-.016}),{.016,.016,.006},0,TILT,.02f,.03f,.06f);dbox(dxf(ct,{0,0,-.0195}),{.01,.01,.001},0,TILT,.15f,.3f,.7f,.5f);
  // АКБ, ремни, разъём, антенна, задний LED
  dbox({0,.0265,.02},{.044,.03,.11},0,0,.14f,.14f,.17f);dbox({0,.0420,.02},{.03,.0008,.06},0,0,.9f,.5f,.1f);
  for(int z=0;z<2;z++)dbox({0,.0265,z?.055:-.015},{.048,.034,.012},0,0,.04f,.04f,.04f);
  dbox({0,.026,.088},{.016,.012,.016},0,0,.9f,.75f,.1f);
  DX at{0,-.5,{0,.035,.12}};dbox(at.c,{.003,.003,.07},0,-.5,.05f,.05f,.05f);dbox(dxf(at,{0,0,.037}),{.006,.006,.006},0,0,1.f,.1f,.1f,.4f);
  float rl=on?1.f:.25f;dbox({0,.012,.087},{.03,.004,.004},0,0,1.f*rl,.1f*rl,.1f*rl,.5f);
  // пропеллеры (прозрачные, рисуются последними)
  glDepthMask(GL_FALSE);
  for(int i=0;i<4;i++){D a=PI/4+i*PI/2,mx=.25*cos(a),mz=.25*sin(a);bool fr=mz<0;D dir=(i%2)?-1:1;
    if(fr)dprop({mx,.03,mz},spin*dir+i*.7,1.f,.5f,.1f);else dprop({mx,.03,mz},spin*dir+i*.7,.7f,.9f,1.f);}
  glDepthMask(GL_TRUE);glLineWidth(1);
  glPopMatrix();}
void renderWorld(V r,V u,V f,V pos,D fov,D asp,bool showDrone,D spin){
  glEnable(GL_DEPTH_TEST);glDepthMask(GL_TRUE);glDisable(GL_TEXTURE_2D);
  glMatrixMode(GL_PROJECTION);glLoadIdentity();persp(fov,asp);glMatrixMode(GL_MODELVIEW);glLoadIdentity();setView(r,u,f,pos);
  glClearColor(.78f,.87f,.95f,1);glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glDisable(GL_FOG);glDepthMask(GL_FALSE);glPushMatrix();glTranslated(pos.x,pos.y,pos.z);sky.draw();
  V sd=norm(Ld)*850,sr=norm(cross(sd,V{0,1,0})),su=norm(cross(sr,sd));glEnable(GL_BLEND);glBlendFunc(GL_SRC_ALPHA,GL_ONE);
  for(int i=3;i>=0;i--){D s=25+i*i*28;glColor4f(1,.93f,.7f,i==0?1.f:.12f);glBegin(GL_QUADS);glVertex3d(sd.x-sr.x*s-su.x*s,sd.y-sr.y*s-su.y*s,sd.z-sr.z*s-su.z*s);glVertex3d(sd.x+sr.x*s-su.x*s,sd.y+sr.y*s-su.y*s,sd.z+sr.z*s-su.z*s);glVertex3d(sd.x+sr.x*s+su.x*s,sd.y+sr.y*s+su.y*s,sd.z+sr.z*s+su.z*s);glVertex3d(sd.x-sr.x*s+su.x*s,sd.y-sr.y*s+su.y*s,sd.z-sr.z*s+su.z*s);glEnd();}
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);glPopMatrix();glDepthMask(GL_TRUE);
  glEnable(GL_FOG);terr.draw();props.draw();
  glColor4f(.16f,.42f,.62f,.82f);glBegin(GL_QUADS);glVertex3d(-400,-3,-400);glVertex3d(400,-3,-400);glVertex3d(400,-3,400);glVertex3d(-400,-3,400);glEnd();
  glDepthMask(GL_FALSE);for(int i=0;i<(int)gates.size();i++)if(i>=gi)drawGate(gates[i],i==gi);
  disc(p.x,hgt(p.x,p.z)+.05,p.z,.6+max(0.,p.y-hgt(p.x,p.z))*.02,.35f);glDepthMask(GL_TRUE);
  if(showDrone)drone(spin);
}

// ---------- 2D: рамка очков, OSD, шумы ----------
void arcPts(D cx,D cy,D r,D a0,D a1,int n,vector<V>&o){for(int i=0;i<=n;i++){D a=a0+(a1-a0)*i/n;o.push_back({cx+r*cos(a),cy+r*sin(a),0});}}
void gogglesFrame(D rx,D ry,D rw,D rh){
  D R=rh*.08;glColor4f(0,0,0,1);
  D cs[4][4]={{rx,ry,PI,1.5*PI},{rx+rw,ry,1.5*PI,2*PI},{rx+rw,ry+rh,0,.5*PI},{rx,ry+rh,.5*PI,PI}};
  D cc[4][2]={{rx+R,ry+R},{rx+rw-R,ry+R},{rx+rw-R,ry+rh-R},{rx+R,ry+rh-R}};
  for(int i=0;i<4;i++){vector<V> pts;arcPts(cc[i][0],cc[i][1],R,cs[i][2],cs[i][3],10,pts);
    glBegin(GL_TRIANGLE_FAN);glVertex2d(cs[i][0],cs[i][1]);for(auto&q2:pts)glVertex2d(q2.x,q2.y);glEnd();}
  vector<V> pts;arcPts(rx+R,ry+R,R,PI,1.5*PI,10,pts);arcPts(rx+rw-R,ry+R,R,1.5*PI,2*PI,10,pts);arcPts(rx+rw-R,ry+rh-R,R,0,.5*PI,10,pts);arcPts(rx+R,ry+rh-R,R,.5*PI,PI,10,pts);
  glLineWidth(3);glColor4f(.35f,.37f,.4f,1);glBegin(GL_LINE_LOOP);for(auto&q2:pts)glVertex2d(q2.x,q2.y);glEnd();
  // видоискатель: угловые скобки
  D m=rh*.045,l=rh*.05;glLineWidth(2);glColor4f(1,1,1,.55f);glBegin(GL_LINES);
  for(int sx=0;sx<2;sx++)for(int sy=0;sy<2;sy++){D x=sx?rx+rw-m:rx+m,y=sy?ry+rh-m:ry+m,dx=sx?-l:l,dy=sy?-l:l;glVertex2d(x,y);glVertex2d(x+dx,y);glVertex2d(x,y);glVertex2d(x,y+dy);}
  glEnd();}
void vignette(D rx,D ry,D rw,D rh){
  D cx=rx+rw/2,cy=ry+rh/2;glBegin(GL_TRIANGLE_STRIP);
  for(int i=0;i<=48;i++){D a=i*2*PI/48;glColor4f(0,0,0,0);glVertex2d(cx+cos(a)*rw*.3,cy+sin(a)*rh*.3);glColor4f(0,0,0,.8f);glVertex2d(cx+cos(a)*rw*.72,cy+sin(a)*rh*.72);}
  glEnd();}
void noiseFx(D rx,D ry,D rw,D rh,D lv){
  int n=(int)(250+lv*7000);glPointSize(2);glBegin(GL_POINTS);
  for(int i=0;i<n;i++){float g=(float)R01();glColor4f(g,g,g,(float)(.2+.6*lv));glVertex2d(rx+R01()*rw,ry+R01()*rh);}glEnd();
  if(lv>.1){glBegin(GL_QUADS);for(int i=0;i<(int)(lv*10);i++){D y=ry+R01()*rh;glColor4f(1,1,1,(float)(.15*lv));rect2(rx,y,rx+rw,y+2+R01()*5);}glEnd();}}
void osd(D rx,D ry,D rw,D rh,V f,V u,V r,D fov){
  D sc=rh/230,m=rh*.085;glPushMatrix();glTranslated(rx,ry,0);
  bool blink=((int)(tm*3))%2==0;D ve=voltage(),rem=1-mah/CAPMAH;bool low=ve<14.0||rem<.2;
  txs(m,m,fmt("%.1fV",ve),sc*1.6,low&&blink?1.f:.78f,low?.3f:1.f,low?.3f:.78f);
  txs(m,m+14*sc,fmt("%dMAH",(int)mah),sc);txs(m,m+24*sc,fmt("%.0fA",curI),sc);
  D hdg=atan2(f.x,-f.z)*180/PI;if(hdg<0)hdg+=360;D cx=rw/2,cy=rh/2;
  glLineWidth(2);glColor4f(.8f,1,.8f,.9f);glBegin(GL_LINES);
  for(int a=-40;a<=40;a++){D d=hdg+a;int dd=((int)floor(d+.5)%360+360)%360;if(dd%10==0){D x=cx+a*rw*.005;glVertex2d(x,m+14*sc);glVertex2d(x,m+(dd%30==0?22:18)*sc);}}
  glEnd();txs(cx-tw(fmt("%03d",(int)hdg),sc*1.3)/2,m-2*sc,fmt("%03d",(int)hdg),sc*1.3);
  glColor4f(1,1,1,1);glBegin(GL_TRIANGLES);glVertex2d(cx,m+24*sc);glVertex2d(cx-4*sc,m+30*sc);glVertex2d(cx+4*sc,m+30*sc);glEnd();
  D dist=hypot(p.x,p.z),lq=clampd(100-95*pow(dist/350,2),0,100);
  txs(rw-m-tw("LQ 100",sc),m,fmt("LQ %d",(int)lq),sc,lq<40?1.f:.78f,lq<40?.4f:1.f,lq<40?.4f:.78f);
  txs(rw-m-tw("SAT 14",sc),m+10*sc,"SAT 14",sc);
  if(((int)(tm*1.5))%2==0){glColor4f(1,.1f,.1f,1);glBegin(GL_TRIANGLE_FAN);for(int i=0;i<12;i++)glVertex2d(rw-m-tw("LQ 100",sc)-8*sc+3*sc*cos(i*PI/6),m+4*sc+3*sc*sin(i*PI/6));glEnd();}
  D spd=hypot(v.x,v.z)*3.6,alt=max(0.,p.y-hgt(p.x,p.z)-.06);
  txs(m,cy-8*sc,fmt("%03d",(int)spd),sc*1.5);txs(m,cy+6*sc,"KMH",sc);
  txs(rw-m-tw("000",sc*1.5),cy-8*sc,fmt("%03d",(int)alt),sc*1.5);txs(rw-m-tw("M",sc),cy+6*sc,"M",sc);
  // горизонт
  D fpx=(rh/2)/tan(fov*PI/360),el=asin(clampd(f.y,-1,1)),roll=atan2(-r.y,u.y)*180/PI;
  glPushMatrix();glTranslated(cx,cy,0);glRotated(-roll,0,0,1);glLineWidth(2);
  for(int a=-60;a<=60;a+=10){D y=fpx*tan(el-a*PI/180);if(fabs(y)>rh*.3)continue;D L=a==0?rw*.2:rw*.07,gap=rw*.04;
    glColor4f(.85f,1,.85f,.85f);glBegin(GL_LINES);glVertex2d(-L-gap,y);glVertex2d(-gap,y);glVertex2d(gap,y);glVertex2d(L+gap,y);
    if(a!=0){D t=a>0?4*sc:-4*sc;glVertex2d(-L-gap,y);glVertex2d(-L-gap,y+t);glVertex2d(L+gap,y);glVertex2d(L+gap,y+t);}glEnd();
    if(a!=0)txs(L+gap+4*sc,y-3*sc,fmt("%d",abs(a)),sc*.8);}
  glPopMatrix();
  glLineWidth(2);glColor4f(1,1,1,.9f);glBegin(GL_LINES);glVertex2d(cx-8*sc,cy);glVertex2d(cx-3*sc,cy);glVertex2d(cx+3*sc,cy);glVertex2d(cx+8*sc,cy);glVertex2d(cx,cy-8*sc);glVertex2d(cx,cy-3*sc);glEnd();
  // низ
  txs(m,rh-m-18*sc,armed?"ARMED":"DISARMED",sc*1.2,armed?.4f:1.f,armed?1.f:.8f,armed?.4f:.3f);
  txs(m,rh-m-6*sc,mode?"ACRO":"ANGLE",sc*1.2);
  D bw=rw*.22,bx=cx-bw/2,by=rh-m-4*sc;glColor4f(1,1,1,.8f);glBegin(GL_LINE_LOOP);glVertex2d(bx,by);glVertex2d(bx+bw,by);glVertex2d(bx+bw,by+8*sc);glVertex2d(bx,by+8*sc);glEnd();
  glColor4f(.5f,1,.5f,.8f);glBegin(GL_QUADS);rect2(bx+2,by+2,bx+2+(bw-4)*thrust,by+8*sc-2);glEnd();
  txs(cx-tw("THR 100",sc)/2,by-10*sc,fmt("THR %d",(int)(thrust*100)),sc);
  int ft=(int)flightT;txs(rw-m-tw("00:00",sc*1.2),rh-m-18*sc,fmt("%02d:%02d",ft/60,ft%60),sc*1.2);
  txs(rw-m-tw("GATE 10/10",sc*1.2),rh-m-6*sc,fmt("GATE %d/%d",min(gi,(int)gates.size()),(int)gates.size()),sc*1.2);
  // стрелка домой
  D bear=atan2(-p.x,p.z)-atan2(f.x,-f.z);glPushMatrix();glTranslated(cx,rh-m-48*sc,0);glRotated(bear*180/PI,0,0,1);
  glColor4f(1,.9f,.2f,.95f);glBegin(GL_TRIANGLES);glVertex2d(0,-9*sc);glVertex2d(-5*sc,6*sc);glVertex2d(5*sc,6*sc);glEnd();glPopMatrix();
  txs(cx+14*sc,rh-m-52*sc,fmt("HOME %dM",(int)dist),sc);
  if(msgT>0&&msg.size())txs(cx-tw(msg,sc*2)/2,rh*.32,msg,sc*2,1,.85f,.3f);
  else if(low&&armed&&blink)txs(cx-tw("LOW BATTERY",sc*2)/2,rh*.32,"LOW BATTERY",sc*2,1,.2f,.2f);
  glPopMatrix();}
void stickBox(D x,D y,D s,D a,D b,const char*l){
  glLineWidth(2);glColor4f(1,1,1,.55f);glBegin(GL_LINE_LOOP);glVertex2d(x,y);glVertex2d(x+s,y);glVertex2d(x+s,y+s);glVertex2d(x,y+s);glEnd();
  glBegin(GL_LINES);glVertex2d(x+s/2,y);glVertex2d(x+s/2,y+s);glVertex2d(x,y+s/2);glVertex2d(x+s,y+s/2);glEnd();
  glColor4f(1,.6f,.1f,1);glBegin(GL_QUADS);rect2(x+s/2+a*s/2-5,y+s/2-b*s/2-5,x+s/2+a*s/2+5,y+s/2-b*s/2+5);glEnd();
  txs(x,y+s+4,l,2,.8f,.8f,.8f);}

// ---------- звук ----------
static void audioCb(void*,Uint8*s,int n){Sint16*o=(Sint16*)s;static D f1=0,f2=0;
  for(int i=0;i<n/2;i++){D f=70+aRpm*300;f1+=f/44100;f2+=f*2.01/44100;
    D vv=((fmod(f1,1)*2-1)*.5+(fmod(f2,1)*2-1)*.25+(R01()-.5)*.25)*aRpm*aRpm*.35;o[i]=(Sint16)(vv*20000);}}

int main(int,char**){
  srand(7);
  if(SDL_Init(SDL_INIT_VIDEO|SDL_INIT_GAMECONTROLLER|SDL_INIT_AUDIO)<0){fprintf(stderr,"SDL: %s\n",SDL_GetError());return 1;}
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE,24);SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS,1);SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES,4);
  win=SDL_CreateWindow("SKYLINE FPV",SDL_WINDOWPOS_CENTERED,SDL_WINDOWPOS_CENTERED,1280,720,SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE);
  if(!win){fprintf(stderr,"Window: %s\n",SDL_GetError());return 1;}
  SDL_GLContext gl=SDL_GL_CreateContext(win);SDL_GL_SetSwapInterval(1);
  SDL_AudioSpec want{};want.freq=44100;want.format=AUDIO_S16SYS;want.channels=1;want.samples=1024;want.callback=audioCb;
  SDL_AudioDeviceID ad=SDL_OpenAudioDevice(nullptr,0,&want,nullptr,0);if(ad)SDL_PauseAudioDevice(ad,0);
  for(int i=0;i<SDL_NumJoysticks();i++)if(SDL_IsGameController(i)){pad=SDL_GameControllerOpen(i);break;}
  buildTerrain();buildProps();
  glEnable(GL_BLEND);glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);glEnable(GL_MULTISAMPLE);
  glFogi(GL_FOG_MODE,GL_LINEAR);GLfloat fc[4]={.78f,.87f,.95f,1};glFogfv(GL_FOG_COLOR,fc);glFogf(GL_FOG_START,160);glFogf(GL_FOG_END,700);
  GLuint tex;glGenTextures(1,&tex);int tw_=0,th_=0;
  printf("SPACE взвод | M режим ANGLE/ACRO | C камера | R заново | P пауза | F полный экран | ESC выход\n");
  Uint64 last=SDL_GetPerformanceCounter();D acc=0,spin=0;bool run=true;V ce={0,5,6};
  while(run){
    SDL_Event e;
    while(SDL_PollEvent(&e)){
      if(e.type==SDL_QUIT)run=false;
      else if(e.type==SDL_KEYDOWN&&!e.key.repeat){switch(e.key.keysym.sym){
        case SDLK_ESCAPE:run=false;break;case SDLK_SPACE:toggleArm();break;case SDLK_m:setMode();break;
        case SDLK_c:camM^=1;break;case SDLK_r:resetAll(true);break;case SDLK_p:paused=!paused;break;
        case SDLK_f:SDL_SetWindowFullscreen(win,(SDL_GetWindowFlags(win)&SDL_WINDOW_FULLSCREEN_DESKTOP)?0:SDL_WINDOW_FULLSCREEN_DESKTOP);break;}}
      else if(e.type==SDL_CONTROLLERDEVICEADDED&&!pad){pad=SDL_GameControllerOpen(e.cdevice.which);}
      else if(e.type==SDL_CONTROLLERDEVICEREMOVED&&pad){SDL_GameControllerClose(pad);pad=nullptr;}
      else if(e.type==SDL_CONTROLLERBUTTONDOWN){switch(e.cbutton.button){
        case SDL_CONTROLLER_BUTTON_A:toggleArm();break;case SDL_CONTROLLER_BUTTON_Y:setMode();break;case SDL_CONTROLLER_BUTTON_X:camM^=1;break;
        case SDL_CONTROLLER_BUTTON_B:resetAll(true);break;case SDL_CONTROLLER_BUTTON_START:paused=!paused;break;}}
    }
    Uint64 now=SDL_GetPerformanceCounter();D dt=min(.05,(now-last)/(D)SDL_GetPerformanceFrequency());last=now;
    gather(dt);
    if(!paused){acc+=dt;while(acc>=1./240){step(1./240);acc-=1./240;}}
    aRpm=(armed&&!paused&&crashT<=0)?.3+.7*thrust:0;spin+=aRpm*60*dt;
    // кадр
    int W,Hh;SDL_GL_GetDrawableSize(win,&W,&Hh);
    glDisable(GL_SCISSOR_TEST);glViewport(0,0,W,Hh);glClearColor(0,0,0,1);glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    bool fpv=camM==0;D rh=fpv?min((D)Hh,W*3./4):Hh,rw=fpv?rh*4/3:W,rx=(W-rw)/2,ry=(Hh-rh)/2;
    D fov=fpv?100:70;V r,u,f,cpos;
    if(fpv){Q jit=axisQ({R01()-.5,R01()-.5,0},.004*thrust*(armed?1:0));Q qc=qnorm(mul(mul(q,axisQ({1,0,0},TILT)),jit));
      r=rot(qc,{1,0,0});u=rot(qc,{0,1,0});f=rot(qc,{0,0,-1});cpos=p+rot(q,{0,.03,-.12});}
    else{V fq=rot(q,{0,0,-1});V fh=norm({fq.x,0,fq.z});V tgt=p-fh*5+V{0,2,0};ce=ce+(tgt-ce)*(1-exp(-dt*6));
      D g=hgt(ce.x,ce.z)+1;if(ce.y<g)ce.y=g;f=norm(p-ce);r=norm(cross(f,V{0,1,0}));u=cross(r,f);cpos=ce;}
    glEnable(GL_SCISSOR_TEST);glScissor((int)rx,(int)ry,(int)rw,(int)rh);glViewport((int)rx,(int)ry,(int)rw,(int)rh);
    renderWorld(r,u,f,cpos,fov,rw/rh,!fpv,spin);
    if(fpv){glBindTexture(GL_TEXTURE_2D,tex);
      if(tw_!=(int)rw||th_!=(int)rh){glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,(int)rw,(int)rh,0,GL_RGB,GL_UNSIGNED_BYTE,nullptr);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);tw_=(int)rw;th_=(int)rh;}
      glCopyTexSubImage2D(GL_TEXTURE_2D,0,0,0,(int)rx,(int)ry,(int)rw,(int)rh);}
    glDisable(GL_SCISSOR_TEST);glViewport(0,0,W,Hh);glDisable(GL_DEPTH_TEST);glDisable(GL_FOG);
    glMatrixMode(GL_PROJECTION);glLoadIdentity();glOrtho(0,W,Hh,0,-1,1);glMatrixMode(GL_MODELVIEW);glLoadIdentity();
    if(fpv){ // бочкообразное искажение линзы FPV-камеры
      glEnable(GL_TEXTURE_2D);glColor3f(1,1,1);const int GX=36,GY=27;const D K=.14;
      for(int j=0;j<GY;j++){glBegin(GL_TRIANGLE_STRIP);
        for(int i=0;i<=GX;i++)for(int s=0;s<2;s++){D sx=i*2./GX-1,sy=(j+s)*2./GY-1,r2=sx*sx+sy*sy,fk=(1+K*r2)/(1+2*K);
          glTexCoord2d(.5+.5*sx*fk,.5+.5*sy*fk);glVertex2d(rx+(sx*.5+.5)*rw,ry+(1-(sy*.5+.5))*rh);}
        glEnd();}
      glDisable(GL_TEXTURE_2D);vignette(rx,ry,rw,rh);}
    D lq=clampd(100-95*pow(hypot(p.x,p.z)/350,2),0,100);
    D nz=clampd((55-lq)/55,0,1);if(crashT>0)nz=1;
    glDisable(GL_SCISSOR_TEST);noiseFx(rx,ry,rw,rh,nz);
    if(fpv)gogglesFrame(rx,ry,rw,rh);
    osd(rx,ry,rw,rh,f,u,r,fov);
    D ss=min(110.,Hh*.15);stickBox(14,Hh-ss-30,ss,in.yaw,in.thr,"YAW/THR");stickBox(W-ss-14,Hh-ss-30,ss,in.rol,in.pit,"ROLL/PIT");
    if(paused)txs(W/2-tw("PAUSED",6)/2,Hh/2,"PAUSED",6,1,1,1);
    SDL_GL_SwapWindow(win);
  }
  if(ad)SDL_CloseAudioDevice(ad);if(pad)SDL_GameControllerClose(pad);
  SDL_GL_DeleteContext(gl);SDL_DestroyWindow(win);SDL_Quit();return 0;}
