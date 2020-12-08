#include<iostream>

int main()
{
    double t=0;
    double dt=0.01;
    int num=0;
    while(t<0.2)
    {
	    num+=1;
	    t+=dt;
	    std::cout<<"This is the "<<num<<"th test : t = "<<t<<" dt = "<<dt<<" t/dt = "<<t/dt<<" int(t/dt) = "<<int(t/dt)<<std::endl;
    }
    //std::stringstream result;
    //result.setf(std::ios::fixed);
    //result.precision(4);
    //result << "u_h_" << int(t/dt) << ".dx";
    //u_h.writeOpenDXData(result.str());

    return 0;
}
