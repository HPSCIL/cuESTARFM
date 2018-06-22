#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "gdal.h"
#include "gdal_priv.h"
#include "gdalwarper.h"
#include <stdio.h>
#include<iostream>
#include"trans.h"
#define num_thread 256
#define num_block 32
__device__ float ftest(int i,int j)
{
	 const float f[]={161,18.51,10.13,7.71,6.61,5.99,5.59,5.32,5.12,4.96,
		4.84,4.75,4.67,4.60,4.54,4.49,4.45,4.41,4.38,4.35,
		4.32,4.30,4.28,4.26,4.24,4.22,4.21,4.20,4.18,4.17,
		4.16,4.15,4.14,4.13,4.121,4.113,4.105,4.098,4.091,4.085,
		4.079,4.073,4.07,4.06,4.06,4.05,4.05,4.043,4.038,4.034
	};
	return f[j];
}
float Stddve(float **a,int n,int width,int height)
{
	float stddve=0,sumx=0,sumxx=0;
	for(int i=0;i<width*height;i++)
	{
		sumx+=a[n][i];
		sumxx+=a[n][i]*a[n][i];
	}
	stddve=sqrt(sumxx/(width*height)-(sumx/(width*height))*(sumx/(width*height)));
	return stddve;
}
__global__ void limit_a_CalcuRela_pairs(float **image_pairs,int num_pairs, int Height,int Width, int Win_size1,float M_err,int BandNum,int current,int *location_p,float *r,float *threshold_d,int task_height  )
{
	const int tid=threadIdx.x;
	const int bid=blockIdx.x;
	const int Idx=num_thread*bid+tid;
	int b=Idx/(Height*Width);
	int j=(Idx-(Height*Width)*b)/Width;
	int i=(Idx-(Height*Width)*b)%Width;
	int k_h=0;
	float dx,dy;
	float sumx,sumy,sumxy,sumxx,sumyy;
	int num=0;
	int ii=0;
	for(int kkk=Idx;kkk<Height*Width;kkk=kkk+num_thread*num_block)
	{
		b=kkk/(Height*Width);
		j=(kkk-(Height*Width)*b)/Width;
		i=(kkk-(Height*Width)*b)%Width;
		dx=0,dy=0;
		num=0;
		sumyy=0;
		sumx=0;
		sumy=0;
		sumxy=0;
		sumxx=0;
		for(int ii=0;ii<BandNum-1;ii++)
		{
			if(image_pairs[ii][j*Width+i]==image_pairs[ii+1][j*Width+i]&&image_pairs[ii+BandNum*num_pairs][j*Width+i]==image_pairs[ii+BandNum*num_pairs+1][j*Width+i])
				num++;
		}

		if(num!=(BandNum-1)||(BandNum==1))
		{
			for( ii=0;ii<BandNum;ii++)
			{
				for(k_h=0;k_h<num_pairs;k_h++)
				{
				sumxy=sumxy+image_pairs[ii+k_h*BandNum][j*Width+i]*image_pairs[ii+(num_pairs+k_h)*BandNum][j*Width+i];
				sumx=sumx+image_pairs[ii+k_h*BandNum][j*Width+i];
				sumy=sumy+image_pairs[ii+(num_pairs+k_h)*BandNum][j*Width+i];
				sumxx=sumxx+image_pairs[ii+k_h*BandNum][j*Width+i]*image_pairs[ii+k_h*BandNum][j*Width+i];
				sumyy=sumyy+image_pairs[ii+(num_pairs+k_h)*BandNum][j*Width+i]*image_pairs[ii+(num_pairs+k_h)*BandNum][j*Width+i];
				}

			}
			dx=sqrt(sumxx/(BandNum*num_pairs)-(sumx/(BandNum*num_pairs))*(sumx/(BandNum*num_pairs)));
			dy=sqrt(sumyy/(BandNum*num_pairs)-(sumy/(BandNum*num_pairs))*(sumy/(BandNum*num_pairs)));
			r[j*Width+i]=(sumxy/(BandNum*num_pairs)-sumx*sumy/(BandNum*BandNum*num_pairs*num_pairs))/(dx*dy);
			if(BandNum==1&&r[j*Width+i]>0)
		      r[j*Width+i]=1;
			if(BandNum==1&&r[j*Width+i]<0)
		      r[j*Width+i]=-1;
		}
		else
		{
			r[j*Width+i]=1;
		}
		if(r[j*Width+i]!=r[j*Width+i])
			r[j*Width+i]=0;
	}
}
__global__ void Blending2_pairs(float **image_pairs,int num_pairs,int Height,int Width, int Win_size1,float M_err,int BandNum,int current,int *location_p,float *r,float *threshold_d,int task_height,float _nodata)
{
	const int tid=threadIdx.x;
	const int bid=blockIdx.x;
	const int Idx=num_thread*bid+tid;
	int b=Idx/(Height*Width);
	int j=(Idx-(Height*Width)*b)/Width;
	int i=(Idx-(Height*Width)*b)%Width;
	int rmin,rmax,smin,smax;
	int r1,s1;
	int Result1=0,m=0;
	int n1;
	float dy;
	float sum,weight_all,weight;
	int k_h=0,k_p=0;
	float aa=0;
	float pix_sum1,pix_sum2;
	//double Aver11;
	//double Aver22;
	float Aver[8];
	float Average1[80],Average3[80];
	float d=0,wi=0;
	float sumx,sumy,sumxy,sumxx,sumyy;
	//float T_weight[8];
	float Aver_all;
	for(int kkk=Idx+current*Width;kkk<(current+task_height)*Width;kkk=kkk+num_thread*num_block)
	{
	
		aa=0;
		pix_sum1=0;
		pix_sum2=0;
		for(m=0;m<80;m++)
		{
			Average1[m]=0;
			Average3[m]=0;
		}
		b=kkk/(Height*Width);
		j=(kkk-(Height*Width)*b)/Width;
		i=(kkk-(Height*Width)*b)%Width;
		for(m=0;m<BandNum;m++)
		{
			pix_sum1+=image_pairs[m][i+Width*j];
			pix_sum2+=image_pairs[m+BandNum][i+Width*j];
		}
		if(fabs(pix_sum1-_nodata)>1e-6&&fabs(pix_sum2-_nodata)>1e-6)
		{
			n1=0;
			weight_all=0,weight=0;
			sum=0;
			sumx=0;
			sumy=0;
		/*	for(m=0;m<8;m++)
			{
				Aver[m]=0;
			}*/
			sumxy=0;
			sumxx=0;
			sumyy=0;
			if(i-Win_size1/2<=0)
				rmin=0;
			else
				rmin = i-Win_size1/2;

			if(i+Win_size1/2>=Width-1)
				rmax = Width-1;
			else
				rmax = i+Win_size1/2;

			if(j-Win_size1/2<=0)
				smin=0;
			else
				smin = j-Win_size1/2;

			if(j+Win_size1/2>=Height-1)
				smax = Height-1;
			else
				smax = j+Win_size1/2;
			r1=rmin,s1=smin;
			for(r1=rmin;r1<=rmax;r1++)
			{
				for(s1=smin;s1<=smax;s1++)
				{  

					Result1=0;	
					for( m=0;m<BandNum*num_pairs;m++)
					{  
						if(fabs(image_pairs[m][ r1+Width*s1]-image_pairs[m][ i+Width*j])<=threshold_d[m])//??
						{
							Result1++;
						}
						else
							break;
					}	

					if(Result1==BandNum*num_pairs )
					{	
						location_p[n1+Idx*Win_size1*Win_size1]=r1+Width*s1;
						d=1+sqrt((float)((r1-i)*(r1-i)+(s1-j)*(s1-j)))/(float)(Win_size1/2);
						weight=1.0/((1.0-r[r1+Width*s1])*d+0.0000001);
						for( m=0;m<BandNum*num_pairs;m++)
						{
<<<<<<< HEAD
							Average1[m]+=(image_pairs[m%BandNum +2*num_pairs*BandNum][r1+Width*s1]-image_pairs[m+num_pairs*BandNum][r1+Width*s1])*weight;
=======
							Average1[m]+=(image_pairs[m%BandNum+2*num_pairs*BandNum][r1+Width*s1]-image_pairs[m+num_pairs*BandNum][r1+Width*s1])*weight;
>>>>>>> 377e39ad22c81bd8790d7be38913158d7e290c00
							//Average2[m]+=(BufferIn55[m][r1+Width*s1]-BufferIn44[m][r1+Width*s1])*weight;
							Average3[m]+=image_pairs[m][r1+Width*s1]*weight;
							//Average4[m]+=BufferIn33[m][r1+Width*s1]*weight;
						}
						weight_all+=weight;
						n1++;
					}
				}
			}

			if(n1>5)
			{
				for(m=0;m<BandNum;m++)
				{
					sumx=0;
					sumy=0;
					sumxy=0;
					sumxx=0;
					sumyy=0;
					for(k_p=0;k_p<8;k_p++)
					{
						Aver[k_p]=0;
					}
						Aver_all=0;
					for(k_h=0;k_h<n1;k_h++)
					{
						for(k_p=0;k_p<num_pairs;k_p++)
						{
						sumxy=sumxy+image_pairs[m+k_p*BandNum][location_p[k_h+Idx*Win_size1*Win_size1]]*image_pairs[m+(k_p+num_pairs)*BandNum][location_p[k_h+Idx*Win_size1*Win_size1]];
						sumx=sumx+image_pairs[m+k_p*BandNum][location_p[k_h+Idx*Win_size1*Win_size1]];
						sumy=sumy+image_pairs[m+(k_p+num_pairs)*BandNum][location_p[k_h+Idx*Win_size1*Win_size1]];
						sumxx=sumxx+image_pairs[m+k_p*BandNum][location_p[k_h+Idx*Win_size1*Win_size1]]*image_pairs[m+k_p*BandNum][location_p[k_h+Idx*Win_size1*Win_size1]];
						sumyy=sumyy+image_pairs[m+(k_p+num_pairs)*BandNum][location_p[k_h+Idx*Win_size1*Win_size1]]*image_pairs[m+(k_p+num_pairs)*BandNum][location_p[k_h+Idx*Win_size1*Win_size1]];
						}
					}
					dy=sqrt(sumyy/(n1*num_pairs)-(sumy/(n1*num_pairs))*(sumy/(n1*num_pairs)));
					if(dy>M_err)
					{
						aa=(sumxy-sumx*sumy/(num_pairs*n1))/(sumyy-sumy*sumy/(n1*num_pairs));
						if(aa>5||aa<0)
						{
							aa=1;
						}
					}
					else
					{
						aa=1.0;
					}
					for( r1=rmin;r1<=rmax;r1++)
					{
						for( s1=smin;s1<=smax;s1++)
						{  
							for(k_p=0;k_p<num_pairs;k_p++)
							{
								Aver[k_p]+=image_pairs[m+2*num_pairs*BandNum][r1+Width*s1]-image_pairs[m+(num_pairs+k_p)*BandNum][r1+Width*s1];
								//Aver22+=BufferIn55[m][r1+Width*s1]-BufferIn44[m][r1+Width*s1];
							}
						}
					}
					for(k_p=0;k_p<num_pairs;k_p++)
					{
						Aver[k_p]=fabs(Aver[k_p])/((float)((rmax-rmin+1)*(smax-smin+1)))+0.0000000001;
						//Aver22+=BufferIn55[m][r1+Width*s1]-BufferIn44[m][r1+Width*s1];
					}
					for(k_p=0;k_p<num_pairs;k_p++)
					{
						Aver_all+=1.0/Aver[k_p];
						//Aver22+=BufferIn55[m][r1+Width*s1]-BufferIn44[m][r1+Width*s1];
					}
					for(k_p=0;k_p<num_pairs;k_p++)
					{
						Aver[k_p]=(1.0/Aver[k_p])/Aver_all;
						//Aver22+=BufferIn55[m][r1+Width*s1]-BufferIn44[m][r1+Width*s1];
					}
					image_pairs[m+(2*num_pairs+1)*BandNum][j*Width+i]=0;
					for(k_p=0;k_p<num_pairs;k_p++)
					{
					image_pairs[m+(2*num_pairs+1)*BandNum][j*Width+i]+=(image_pairs[m+k_p*BandNum][j*Width+i]+aa*Average1[m+k_p*BandNum]/weight_all)*Aver[k_p];
					}
					if(image_pairs[m+(2*num_pairs+1)*BandNum][j*Width+i]<0||image_pairs[m+(2*num_pairs+1)*BandNum][j*Width+i]>1)
					{
						image_pairs[m+(2*num_pairs+1)*BandNum][j*Width+i]=0;
						for(k_p=0;k_p<num_pairs;k_p++)
						{
							image_pairs[m+(2*num_pairs+1)*BandNum][j*Width+i]+=Average3[m+k_p*BandNum]*Aver[k_p]/weight_all;
						}
					}
				}
			}
			else
			{
				for(m=0;m<BandNum;m++)
				{
					for(k_p=0;k_p<8;k_p++)
					{
						Aver[k_p]=0;
					}
					Aver_all=0;
					for( r1=rmin;r1<=rmax;r1++)
					{
						for( s1=smin;s1<=smax;s1++)
						{  
							for(k_p=0;k_p<num_pairs;k_p++)
							{
								Aver[k_p]+=image_pairs[m+2*num_pairs*BandNum][r1+Width*s1]-image_pairs[m+(num_pairs+k_p)*BandNum][r1+Width*s1];
								//Aver22+=BufferIn55[m][r1+Width*s1]-BufferIn44[m][r1+Width*s1];
							}
						}
					}
					for(k_p=0;k_p<num_pairs;k_p++)
					{
						Aver[k_p]=fabs(Aver[k_p])/((float)((rmax-rmin+1)*(smax-smin+1)))+0.0000000001;
						//Aver22+=BufferIn55[m][r1+Width*s1]-BufferIn44[m][r1+Width*s1];
					}
					for(k_p=0;k_p<num_pairs;k_p++)
					{
						Aver_all+=1.0/Aver[k_p];
						//Aver22+=BufferIn55[m][r1+Width*s1]-BufferIn44[m][r1+Width*s1];
					}
					for(k_p=0;k_p<num_pairs;k_p++)
					{
						Aver[k_p]=(1.0/Aver[k_p])/Aver_all;
						//Aver22+=BufferIn55[m][r1+Width*s1]-BufferIn44[m][r1+Width*s1];
					}	
					image_pairs[m+(2*num_pairs+1)*BandNum][j*Width+i]=0;
					for(k_p=0;k_p<num_pairs;k_p++)
					{
						image_pairs[m+(2*num_pairs+1)*BandNum][j*Width+i]+=image_pairs[m+k_p*BandNum][j*Width+i]*Aver[k_p];
					}
				}
			}
		}
		else
		{
			for(m=0;m<BandNum;m++)
			{
				image_pairs[m+(2*num_pairs+1)*BandNum][j*Width+i]=0;
			}
		}
	}
}
__global__ void limit_a_CalcuRela(float **BufferIn11,float **BufferIn22,float **BufferIn33,float **BufferIn44,float **BufferIn55, int Height,int Width, int Win_size1,float M_err,int BandNum,int current,int *location_p,float *r,float *threshold_d,int task_height,float _nodata  )
{
	const int tid=threadIdx.x;
	const int bid=blockIdx.x;
	const int Idx=num_thread*bid+tid;
	int b=Idx/(Height*Width);
	int j=(Idx-(Height*Width)*b)/Width;
	int i=(Idx-(Height*Width)*b)%Width;
	int r1,s1;
	int k_h=0;
	float dx,dy;
	float sumx,sumy,sumxy,sumxx,sumyy;
	int num=0;
	for(int kkk=Idx;kkk<Height*Width;kkk=kkk+num_thread*num_block)
	{
		b=kkk/(Height*Width);
		j=(kkk-(Height*Width)*b)/Width;
		i=(kkk-(Height*Width)*b)%Width;
		dx=0,dy=0;
		num=0;
		sumyy=0;
		sumx=0;
		sumy=0;
		sumxy=0;
		sumxx=0;
		for(int ii=0;ii<BandNum-1;ii++)
		{
			if(BufferIn11[ii][j*Width+i]==BufferIn11[ii+1][j*Width+i]&&BufferIn22[ii][j*Width+i]==BufferIn22[ii+1][j*Width+i])
				num++;
		}

		if(num!=(BandNum-1)||(BandNum==1))
		{
			for(int ii=0;ii<BandNum;ii++)
			{
				sumxy=sumxy+BufferIn11[ii][j*Width+i]*BufferIn22[ii][j*Width+i]+BufferIn33[ii][j*Width+i]*BufferIn44[ii][j*Width+i];
				sumx=sumx+BufferIn11[ii][j*Width+i]+BufferIn33[ii][j*Width+i];
				sumy=sumy+BufferIn22[ii][j*Width+i]+BufferIn44[ii][j*Width+i];
				sumxx=sumxx+BufferIn11[ii][j*Width+i]*BufferIn11[ii][j*Width+i]+BufferIn33[ii][j*Width+i]*BufferIn33[ii][j*Width+i];
				sumyy=sumyy+BufferIn22[ii][j*Width+i]*BufferIn22[ii][j*Width+i]+BufferIn44[ii][j*Width+i]*BufferIn44[ii][j*Width+i];

			}
			dx=sqrt(sumxx/(BandNum*2)-(sumx/(BandNum*2))*(sumx/(BandNum*2)));
			dy=sqrt(sumyy/(BandNum*2)-(sumy/(BandNum*2))*(sumy/(BandNum*2)));
			r[j*Width+i]=(sumxy/(BandNum*2)-sumx*sumy/(BandNum*BandNum*4))/(dx*dy);
			if(BandNum==1&&r[j*Width+i]>0)
		      r[j*Width+i]=1;
			if(BandNum==1&&r[j*Width+i]<0)
		      r[j*Width+i]=-1;
		}
		else
		{
			r[j*Width+i]=1;
		}
		if(r[j*Width+i]!=r[j*Width+i])
			r[j*Width+i]=0;
	}

	
}
__global__ void Blending2(float **BufferIn11,float **BufferIn22,float **BufferIn33,float **BufferIn44,float **BufferIn55, float **BufferOut,int Height,int Width, int Win_size1,float M_err,int BandNum,int current,int *location_p,float *r,float *threshold_d,int task_height ,float _nodata)
{
	const int tid=threadIdx.x;
	const int bid=blockIdx.x;
	const int Idx=num_thread*bid+tid;
	int b=Idx/(Height*Width);
	int j=(Idx-(Height*Width)*b)/Width;
	int i=(Idx-(Height*Width)*b)%Width;
	int rmin,rmax,smin,smax;
	int r1,s1;
	int Result1=0,m=0;
	int n1;
	float dy;
	float sum,weight_all,weight;
	int k_h=0;
	float aa=0;
	float pix_sum1,pix_sum2;
	double Aver11;
	double Aver22;
	float Average1[10],Average2[10],Average3[10],Average4[10];
	float d=0,wi=0;
	float sumx,sumy,sumxy,sumxx,sumyy;
	float T1_weight;
	float T2_weight;
	for(int kkk=Idx+current*Width;kkk<(current+task_height)*Width;kkk=kkk+num_thread*num_block)
	{
		aa=0;
		pix_sum1=0;
		pix_sum2=0;
		for(m=0;m<10;m++)
		{
			Average1[m]=0;
			Average2[m]=0;
			Average3[m]=0;
			Average4[m]=0;
		}
		b=kkk/(Height*Width);
		j=(kkk-(Height*Width)*b)/Width;
		i=(kkk-(Height*Width)*b)%Width;
		for(m=0;m<BandNum;m++)
		{
			pix_sum1+=BufferIn11[m][i+Width*j];
			pix_sum2+=BufferIn33[m][i+Width*j];
		}
		if(fabs(pix_sum1-_nodata)>1e-6&&fabs(pix_sum2-_nodata)>1e-6) 
		{
			n1=0;
			weight_all=0,weight=0;
			sum=0;
			sumx=0;
			Aver11=0;
			Aver22=0;
			sumy=0;
			sumxy=0;
			sumxx=0;
			sumyy=0;
			if(i-Win_size1/2<=0)
				rmin=0;
			else
				rmin = i-Win_size1/2;

			if(i+Win_size1/2>=Width-1)
				rmax = Width-1;
			else
				rmax = i+Win_size1/2;

			if(j-Win_size1/2<=0)
				smin=0;
			else
				smin = j-Win_size1/2;

			if(j+Win_size1/2>=Height-1)
				smax = Height-1;
			else
				smax = j+Win_size1/2;
			r1=rmin,s1=smin;
			for(r1=rmin;r1<=rmax;r1++)
			{
				for(s1=smin;s1<=smax;s1++)
				{  

					Result1=0;	
					for( m=0;m<BandNum;m++)
					{  
						if(fabs(BufferIn11[m][ r1+Width*s1]-BufferIn11[m][ i+Width*j])<=threshold_d[m]&&fabs(BufferIn33[m][ r1+Width*s1]-BufferIn33[m][ i+Width*j])<=threshold_d[m+BandNum])//??
						{
							Result1++;
						}
						else
							break;
					}	

					if(Result1==BandNum )
					{	
						location_p[n1+Idx*Win_size1*Win_size1]=r1+Width*s1;
						d=1+sqrt((float)((r1-i)*(r1-i)+(s1-j)*(s1-j)))/(float)(Win_size1/2);
<<<<<<< HEAD
						weight=1.0/((1.0-r[r1+Width*s1])*d+0.0000001 );
=======
						weight=1.0/((1.0-r[r1+Width*s1])*d+0.0000001);
>>>>>>> 377e39ad22c81bd8790d7be38913158d7e290c00
						for( m=0;m<BandNum;m++)
						{
							Average1[m]+=(BufferIn55[m][r1+Width*s1]-BufferIn22[m][r1+Width*s1])*weight;
							Average2[m]+=(BufferIn55[m][r1+Width*s1]-BufferIn44[m][r1+Width*s1])*weight;
							Average3[m]+=BufferIn11[m][r1+Width*s1]*weight;
							Average4[m]+=BufferIn33[m][r1+Width*s1]*weight;
						}
						weight_all+=weight;
						n1++;
					}
				}
			}

			if(n1>5)
			{
				for(m=0;m<BandNum;m++)
				{
					sumx=0;
					sumy=0;
					sumxy=0;
					sumxx=0;
					sumyy=0;
					Aver11=0;
					Aver22=0;
					for(k_h=0;k_h<n1;k_h++)
					{
						sumxy=sumxy+BufferIn11[m][location_p[k_h+Idx*Win_size1*Win_size1]]*BufferIn22[m][location_p[k_h+Idx*Win_size1*Win_size1]]+BufferIn33[m][location_p[k_h+Idx*Win_size1*Win_size1]]*BufferIn44[m][location_p[k_h+Idx*Win_size1*Win_size1]];
						sumx=sumx+BufferIn11[m][location_p[k_h+Idx*Win_size1*Win_size1]]+BufferIn33[m][location_p[k_h+Idx*Win_size1*Win_size1]];
						sumy=sumy+BufferIn22[m][location_p[k_h+Idx*Win_size1*Win_size1]]+BufferIn44[m][location_p[k_h+Idx*Win_size1*Win_size1]];
						sumxx=sumxx+BufferIn11[m][location_p[k_h+Idx*Win_size1*Win_size1]]*BufferIn11[m][location_p[k_h+Idx*Win_size1*Win_size1]]+BufferIn33[m][location_p[k_h+Idx*Win_size1*Win_size1]]*BufferIn33[m][location_p[k_h+Idx*Win_size1*Win_size1]];
						sumyy=sumyy+BufferIn22[m][location_p[k_h+Idx*Win_size1*Win_size1]]*BufferIn22[m][location_p[k_h+Idx*Win_size1*Win_size1]]+BufferIn44[m][location_p[k_h+Idx*Win_size1*Win_size1]]*BufferIn44[m][location_p[k_h+Idx*Win_size1*Win_size1]];
					}
					dy=sqrt(sumyy/(n1*2)-(sumy/(n1*2))*(sumy/(n1*2)));
					if(dy>M_err)
					{
						aa=(sumxy-sumx*sumy/(2*n1))/(sumyy-sumy*sumy/(n1*2));
						if(aa>5||aa<0)
						{
							aa=1;
						}
					}
					else
					{
						aa=1.0;
					}
					for(int r1=rmin;r1<=rmax;r1++)
					{
						for(int s1=smin;s1<=smax;s1++)
						{  
							Aver11+=BufferIn55[m][r1+Width*s1]-BufferIn22[m][r1+Width*s1];
							Aver22+=BufferIn55[m][r1+Width*s1]-BufferIn44[m][r1+Width*s1];	
						}
					}
					Aver11=fabs(Aver11)/((float)((rmax-rmin+1)*(smax-smin+1)))+0.0000000001;
					Aver22=fabs(Aver22)/((float)((rmax-rmin+1)*(smax-smin+1)))+0.0000000001;
					T1_weight=1.0/Aver11/(1.0/Aver11+1.0/Aver22);
					T2_weight=1.0/Aver22/(1.0/Aver11+1.0/Aver22);	
					BufferOut[m][j*Width+i]=(BufferIn11[m][j*Width+i]+aa*Average1[m]/weight_all)*T1_weight+(BufferIn33[m][j*Width+i]+aa*Average2[m]/weight_all)*T2_weight;
					if(BufferOut[m][j*Width+i]<0)
						BufferOut[m][j*Width+i]=Average3[m]*T1_weight/weight_all+Average4[m]*T2_weight/weight_all;
				}
			}
			else
			{
				for(m=0;m<BandNum;m++)
				{
					Aver11=0;
					Aver22=0;
					for(int r1=rmin;r1<=rmax;r1++)
					{
						for(int s1=smin;s1<=smax;s1++)
						{  
							Aver11+=BufferIn55[m][r1+Width*s1]-BufferIn22[m][r1+Width*s1];
							Aver22+=BufferIn55[m][r1+Width*s1]-BufferIn44[m][r1+Width*s1];		
						}
					}
					Aver11=fabs(Aver11)/((float)((rmax-rmin+1)*(smax-smin+1)))+0.0000000001;
					Aver22=fabs(Aver22)/((float)((rmax-rmin+1)*(smax-smin+1)))+0.0000000001;
					T1_weight=1/Aver11/(1/Aver11+1/Aver22);
					T2_weight=1/Aver22/(1/Aver11+1/Aver22);	
					BufferOut[m][j*Width+i]=BufferIn11[m][j*Width+i]*T1_weight+BufferIn33[m][j*Width+i]*T2_weight;
				}
			}
		}
		else
		{
			for(m=0;m<BandNum;m++)
			{
				BufferOut[m][j*Width+i]=0;
			}
		}
	}
}
void runtest1_pairs(float **sub_area,int Height,int Width,PARAMETER *par,int BandNum,float *std,int current,int task_height)
 {
	 float **dev_sub_area;
	 float **a;
	 float *dev_std,*r;
	 int num_pairs=par->NUM_PAIRS;
	 int windows=par->WIN_SIZE;
	 float _nodata=par->_nodata;
	 float M_err=par->uncertain;
	 a = (float**)malloc(2*(par->NUM_PAIRS+1)*BandNum*sizeof(float*));
	 for(int b=0;b<2*(par->NUM_PAIRS+1)*BandNum;b++)
	 {
		 cudaMalloc((void**)&a[b],Height*Width*sizeof(float));
	 }
	 int *Location_P;
	cudaMalloc((void***)&dev_sub_area,sizeof(float*)*2*(par->NUM_PAIRS+1)*BandNum);
	cudaMalloc((void**)&Location_P,sizeof(float)*par->WIN_SIZE*par->WIN_SIZE*num_block*num_thread);
	cudaMalloc((void**)&r,sizeof(float)*Height*Width);
	cudaMalloc((void**)&dev_std,sizeof(float)*BandNum*par->NUM_PAIRS);
	 for(int b=0;b<2*(par->NUM_PAIRS+1)*BandNum;b++)
	 {
        cudaMemcpy(a[b], sub_area[b],Height*Width*sizeof(float),cudaMemcpyHostToDevice);
	 }
	 cudaMemcpy(dev_std, std,sizeof(float)*BandNum*par->NUM_PAIRS,cudaMemcpyHostToDevice);
	 cudaMemcpy(dev_sub_area,a,sizeof(float*)*2*(par->NUM_PAIRS+1)*BandNum,cudaMemcpyHostToDevice);
	 limit_a_CalcuRela_pairs<<<num_block, num_thread>>>(dev_sub_area, num_pairs,Height, Width,windows,M_err,BandNum,current,Location_P,r,dev_std,task_height);
	// cudaMemGetInfo(&ff, &tt);
	Blending2_pairs<<<num_block, num_thread>>>(dev_sub_area,num_pairs, Height, Width,  windows, M_err,BandNum,current,Location_P,r,dev_std,task_height,_nodata);
	for(int g=0;g<BandNum;g++)
	{
		cudaMemcpy(sub_area[g+(2*num_pairs+1)*BandNum],a[g+(2*num_pairs+1)*BandNum],Height*Width*sizeof(float),cudaMemcpyDeviceToHost);
	}
	for(int g=0;g<BandNum;g++)
	{
		cudaFree(a[g]);
	}
	cudaFree(Location_P);
	cudaFree(r);
	cudaFree(dev_std);
 }
void runtest1(float **BufferIn11,float **BufferIn22,float **BufferIn33,float **BufferIn44,float **BufferIn55,float **BufferOut,int Height,int Width,int Win_size1,float M_err,int BandNum,float *std,int current,int task_height,float _nodata)
{
	float **dev_BufferIn11,**dev_BufferIn22,**dev_BufferIn33,**dev_BufferIn44,**dev_BufferIn55,**dev_BufferOut;//*Changed_BufferIn11,*Changed_BufferIn33;
	float **a,**f,**c,**d,**e,**out;
	//unsigned int ff, tt
	float *dev_std,*r;
	a = (float**)malloc(BandNum*sizeof(float*));
	f = (float**)malloc(BandNum*sizeof(float*));
	c = (float**)malloc(BandNum*sizeof(float*));
	d = (float**)malloc(BandNum*sizeof(float*));
	e = (float**)malloc(BandNum*sizeof(float*));
	out=(float**)malloc(BandNum*sizeof(float*));
	//test=(float*)malloc(BandNum*Height*Width*sizeof(float));
	for(int b=0;b<BandNum;b++)
	{
		cudaMalloc((void**)&a[b],Height*Width*sizeof(float));
		cudaMalloc((void**)&f[b],Height*Width*sizeof(float));
		cudaMalloc((void**)&c[b],Height*Width*sizeof(float));
		cudaMalloc((void**)&d[b],Height*Width*sizeof(float));
		cudaMalloc((void**)&e[b],Height*Width*sizeof(float));
		cudaMalloc((void**)&out[b],Height*Width*sizeof(float));

	}
	//int num_block= Height* Width*BandNum/num_thread+1;
	int *Location_P;
	cudaMalloc((void***)&dev_BufferIn11,sizeof(float*)*BandNum);
	cudaMalloc((void***)&dev_BufferIn22,sizeof(float*)*BandNum);
	cudaMalloc((void***)&dev_BufferIn33,sizeof(float*)*BandNum);
	cudaMalloc((void***)&dev_BufferIn44,sizeof(float*)*BandNum);
	cudaMalloc((void***)&dev_BufferIn55,sizeof(float*)*BandNum);
	cudaMalloc((void***)&dev_BufferOut,sizeof(float*)*BandNum);
	cudaMalloc((void**)&Location_P,sizeof(float)*Win_size1*Win_size1*num_block*num_thread);
	//cudaMalloc((void**)&Changed_BufferIn11,sizeof(float)*Height*Width*BandNum);
	//cudaMalloc((void**)&Changed_BufferIn33,sizeof(float)*Height*Width*BandNum);
	cudaMalloc((void**)&r,sizeof(float)*Height*Width);
	cudaMalloc((void**)&dev_std,sizeof(float)*BandNum*2);
	//cudaMalloc((void**)&Location_P,sizeof(float)*100*Width*Height*BandNum);
	for(int g=0;g<BandNum;g++)
	{
		cudaMemcpy(a[g], BufferIn11[g],Height*Width*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(f[g], BufferIn22[g],Height*Width*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(c[g], BufferIn33[g],Height*Width*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(d[g], BufferIn44[g],Height*Width*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(e[g], BufferIn55[g],Height*Width*sizeof(float),cudaMemcpyHostToDevice);
	}
	cudaMemcpy(dev_BufferIn11,a,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_BufferIn22,f,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_BufferIn33, c,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_BufferIn44, d,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_BufferIn55, e,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_BufferOut,out,sizeof(float*)*BandNum,cudaMemcpyHostToDevice);
	cudaMemcpy(dev_std, std,sizeof(float)*BandNum*2,cudaMemcpyHostToDevice);
	int deviceCount;
	cudaGetDeviceCount(&deviceCount);
	//cudaSetDevice(0);
	/*size_t ff,tt;
    cudaMemGetInfo(&ff, &tt);
	cudaEvent_t start,stop;
	cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start,0);*/
	
	limit_a_CalcuRela<<<num_block, num_thread>>>(dev_BufferIn11,dev_BufferIn22,dev_BufferIn33,dev_BufferIn44,dev_BufferIn55, Height, Width,  Win_size1,M_err,BandNum,current,Location_P,r,dev_std,task_height,_nodata);
	// cudaMemGetInfo(&ff, &tt);
	Blending2<<<num_block, num_thread>>>(dev_BufferIn11,dev_BufferIn22,dev_BufferIn33,dev_BufferIn44,dev_BufferIn55,dev_BufferOut, Height, Width,  Win_size1, M_err,BandNum,current,Location_P,r,dev_std,task_height,_nodata);
	/*cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
    float costtime=0.0f;
    cudaEventElapsedTime(&costtime,start,stop);
	std::cout<<costtime/1000<<"  ";*/
//	 cudaMemGetInfo(&ff, &tt);
	for(int g=0;g<BandNum;g++)
	{
		cudaMemcpy(BufferOut[g],out[g],Height*Width*sizeof(float),cudaMemcpyDeviceToHost);
	}
	//cudaMemcpy(test,Changed_BufferIn11,sizeof(float)*BandNum*Height*Width,cudaMemcpyDeviceToHost);
	//for(int i=0;i<100;i++)
	//std::cout<<test[i]<<" ";
	/*if(BufferOut[1][1]<0)
		std::cout<<"wrong";*/
	for(int g=0;g<BandNum;g++)
	{
		cudaFree(a[g]);
		cudaFree(f[g]);
		cudaFree(c[g]);
		cudaFree(d[g]);
		cudaFree(e[g]);
		cudaFree(out[g]);
	}
	//cudaFree(Changed_BufferIn11);
	//cudaFree(Changed_BufferIn33);
	cudaFree(Location_P);
	cudaFree(r);
	cudaFree(dev_std);
}
void runtest(float **BufferIn11,float **BufferIn22,float **BufferIn33,float **BufferIn44,float **BufferIn55,float **BufferOut,int Height,int Width,int Win_size1,float M_err,int num_class,int BandNum,float _nodata)
{
	float *std=new float[BandNum*2];
	for(int i=0;i<BandNum;i++)
	{
		std[i]=Stddve(BufferIn11,i,Width,Height)*2.0/num_class;
		std[i+BandNum]=Stddve(BufferIn33,i,Width,Height)*2.0/num_class;
		std::cout<<std[i]<<"  "<<std[i+BandNum]<<"  ";
	}
	int maxnum;
	size_t ff,tt;
	//cudaSetDevice(0);
	cudaMemGetInfo(&ff, &tt);
	maxnum=(ff-sizeof(float)*Win_size1*Win_size1*num_block*num_thread*2)/(BandNum*sizeof(float)*8);
	int sub_height=maxnum/Width-Win_size1;
	int kk=0;
	int i,j;
	float **sub_BufferIn11,**sub_BufferIn22,**sub_BufferIn33,**sub_BufferIn44,**sub_BufferIn55,**sub_out;
	for(int heiht_all=0;heiht_all<Height;heiht_all+=sub_height)
	{
		int task_start=kk*sub_height;
		int task_end;
		if((kk+1)*sub_height-Height<=0)
			task_end=(kk+1)*sub_height-1;
		else
			task_end=Height-1; 
		int data_start,data_end;
		if(task_start-Win_size1/2<=0)
			data_start= 0;
		else
			data_start=task_start-Win_size1/2;
		if(task_end+Win_size1/2>=Height-1)
			data_end=Height-1;
		else
			data_end=task_end+Win_size1/2;
		int data_height=data_end-data_start+1;
		sub_BufferIn11=(float**)malloc(BandNum*sizeof(float*));
		sub_BufferIn22=(float**)malloc(BandNum*sizeof(float*));
		sub_BufferIn33=(float**)malloc(BandNum*sizeof(float*));
		sub_BufferIn44=(float**)malloc(BandNum*sizeof(float*));
		sub_BufferIn55=(float**)malloc(BandNum*sizeof(float*));
		sub_out=(float**)malloc(BandNum*sizeof(float*));
		for(int b=0;b<BandNum;b++)
		{
			sub_BufferIn11[b]=new float[data_height*Width];
			sub_BufferIn22[b]=new float[data_height*Width];
			sub_BufferIn33[b]=new float[data_height*Width];
			sub_BufferIn44[b]=new float[data_height*Width];
			sub_BufferIn55[b]=new float[data_height*Width];
			sub_out[b]=new float[data_height*Width];
		}
		int copy;
		for(int k=0;k<BandNum;k++)
		{
			copy=0;
			for( i=data_start;i<=data_end;i++)
			{
				for( j=0;j<Width;j++)
				{
					sub_BufferIn11[k][copy*Width+j]=BufferIn11[k][i*Width+j];
					sub_BufferIn22[k][copy*Width+j]=BufferIn22[k][i*Width+j];
					sub_BufferIn33[k][copy*Width+j]=BufferIn33[k][i*Width+j];
					sub_BufferIn44[k][copy*Width+j]=BufferIn44[k][i*Width+j];
					sub_BufferIn55[k][copy*Width+j]=BufferIn55[k][i*Width+j];
				}
				copy++;
			}
		}
		int current=task_start-data_start;
		int task_height=task_end-task_start+1;
		runtest1(sub_BufferIn11,sub_BufferIn22,sub_BufferIn33,sub_BufferIn44,sub_BufferIn55,sub_out,data_height,Width,Win_size1,M_err, BandNum,std,current,task_height,_nodata);
		
		for(int k=0;k<BandNum;k++)
		{
			current=task_start-data_start;
			for(int i=task_start;i<=task_end;i++)
			{
				for(int j=0;j<Width;j++)
				{
					BufferOut[k][i*Width+j]=sub_out[k][current*Width+j];
				}
				current++;
			}
		}
		for(int g=0;g<BandNum;g++)
	{
		delete sub_BufferIn11[g];
		delete sub_BufferIn22[g];
		delete sub_BufferIn33[g];
		delete sub_BufferIn44[g];
		delete sub_BufferIn55[g];
		delete sub_out[g];
		/*cudaFree(dev_BufferIn11[g]);
		cudaFree(dev_BufferIn22[g]);
		cudaFree(dev_BufferIn33[g]);
		cudaFree(dev_BufferIn44[g]);
		cudaFree(dev_BufferIn55[g]);
		cudaFree(dev_BufferOut[g]);*/
	}
		kk++;
	}
}
void runtest_pairs(CuLayer *psensor,PARAMETER *par,int solve)
{
	int Height=psensor[0].getHeight();
	int Width=psensor[0].getWidth();
	int BandNum=psensor[0].getbandCount();
	float *std=new float[psensor[0].getbandCount()*par->NUM_PAIRS];
	for(int j=0;j<par->NUM_PAIRS;j++)
	{
		for(int i=0;i<psensor[0].getbandCount();i++)
		{

			std[i+psensor[0].getbandCount()*j]=Stddve(psensor[j].getData(),i,psensor[0].getWidth(),psensor[0].getHeight())*2.0/par->class_num;
			std::cout<<std[i+psensor[0].getbandCount()*j]<<"  ";
		}
	}
	int maxnum;
	size_t ff,tt;
	//cudaSetDevice(0);
	cudaMemGetInfo(&ff, &tt);
	maxnum=(ff-sizeof(float)*par->WIN_SIZE*par->WIN_SIZE*num_block*num_thread*2)/(psensor[0].getbandCount()*sizeof(float)*2*(par->NUM_PAIRS+1))-par->WIN_SIZE;
	int sub_height=maxnum/psensor[0].getWidth();
	int kk=0;
	int i,j,c;
	float **sub_area;
	for(int heiht_all=0;heiht_all<Height;heiht_all+=sub_height)
	{
		int task_start=kk*sub_height;
		int task_end;
		if((kk+1)*sub_height-Height<=0)
			task_end=(kk+1)*sub_height-1;
		else
			task_end=Height-1; 
		int data_start,data_end;
		if(task_start-par->WIN_SIZE/2<=0)
			data_start= 0;
		else
			data_start=task_start-par->WIN_SIZE/2;
		if(task_end+par->WIN_SIZE/2>=Height-1)
			data_end=Height-1;
		else
			data_end=task_end+par->WIN_SIZE/2;
		int data_height=data_end-data_start+1;
		sub_area=(float**)malloc(2*(par->NUM_PAIRS+1)*BandNum*sizeof(float*));
		for(int b=0;b<2*(par->NUM_PAIRS+1)*BandNum;b++)
		{
			sub_area[b]=new float[data_height*Width];
		}
		int copy;
		for(int k=0;k<BandNum;k++)
		{
			copy=0;
			for( i=data_start;i<=data_end;i++)
			{
				for( j=0;j<Width;j++)
				{
					for(c=0;c<2*par->NUM_PAIRS;c++)
					{
						sub_area[k+c*BandNum][copy*Width+j]=psensor[c].getData()[k][i*Width+j];
					}
					sub_area[2*par->NUM_PAIRS*BandNum+k][copy*Width+j]=psensor[2*(par->NUM_PAIRS+solve)].getData()[k][i*Width+j];
				}
				copy++;
			}
		}
		int current=task_start-data_start;
		int task_height=task_end-task_start+1;
		runtest1_pairs(sub_area,data_height,Width,par, BandNum,std,current,task_height);
		
		for(int k=0;k<BandNum;k++)
		{
			current=task_start-data_start;
			for(int i=task_start;i<=task_end;i++)
			{
				for(int j=0;j<Width;j++)
				{
					psensor[2*(par->NUM_PAIRS+solve)+1].getData()[k][i*Width+j]=sub_area[(2*par->NUM_PAIRS+1)*BandNum+k][current*Width+j];
				}
				current++;
			}
		}
		for(int b=0;b<2*(par->NUM_PAIRS+1)*BandNum;b++)
		{
			delete sub_area[b];
			/*cudaFree(dev_BufferIn11[g]);
			cudaFree(dev_BufferIn22[g]);
			cudaFree(dev_BufferIn33[g]);
			cudaFree(dev_BufferIn44[g]);
			cudaFree(dev_BufferIn55[g]);
			cudaFree(dev_BufferOut[g]);*/
		}
		kk++;
	}
}
void Re_fusion3(CuLayer *psensor,PARAMETER *par)
{
	int i,j,m,c;
	long now1 = clock();
	for(c=0;c<par->NUM_PREDICTIONS;c++)
	{
		psensor[2*(par->NUM_PAIRS+c)+1].resize(psensor[0].getWidth(),psensor[0].getHeight(),psensor[0].getbandCount());
		if(par->NUM_PAIRS==2)
		{
			runtest(psensor[0].getData(),psensor[2].getData(),psensor[1].getData(),psensor[3].getData(),psensor[2*(par->NUM_PAIRS+c)].getData(),psensor[2*(par->NUM_PAIRS+c)+1].getData(),psensor[0].getHeight(),psensor[0].getWidth(),par->WIN_SIZE,par->uncertain,par->class_num,psensor[0].getbandCount(),par->_nodata);
			//char* driverName = "GTiff";
			psensor[2*(par->NUM_PAIRS+c)+1].setGeoTransform(psensor[0].getGeoTransform());
			psensor[2*(par->NUM_PAIRS+c)+1].setProjection(psensor[0].getProjection());
			psensor[2*(par->NUM_PAIRS+c)+1].Write(psensor[2*(par->NUM_PAIRS+c)+1].outpath,par->G_Type);
		}
		else
		{
			runtest_pairs(psensor,par,c);
			psensor[2*(par->NUM_PAIRS+c)+1].setGeoTransform(psensor[0].getGeoTransform());
			psensor[2*(par->NUM_PAIRS+c)+1].setProjection(psensor[0].getProjection());
			psensor[2*(par->NUM_PAIRS+c)+1].Write(psensor[2*(par->NUM_PAIRS+c)+1].outpath,par->G_Type);
		}
	}

}
 void Re_fusion2(const char * BufferIn0,const char * BufferIn1,const char * BufferIn2,const char * BufferIn3,const char * BufferIn4,const char * BufferOut,int win_size,int class_num,float M_err,float _nodata)
{
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO"); 
	GDALDataset *Landsat0 = (GDALDataset*) GDALOpen(BufferIn0,GA_ReadOnly);
	int width,height,BandNum;
	width = Landsat0->GetRasterXSize();
	height = Landsat0->GetRasterYSize();
	BandNum = Landsat0->GetRasterCount();
	//height=2000;
   // width=2000;
	float** BufferLandsat_0 = new float*[BandNum];
	int b,k;
	for( b=0;b<BandNum;b++)
	{
		BufferLandsat_0[b] = new float[width*height];
	//	printf("%u\n", BufferLandsat_0[b]);
	}
	
	for( k=0;k<BandNum;k++)
	{
		GDALRasterBand* hInBand1 = Landsat0->GetRasterBand(k+1);
		hInBand1->RasterIO(GF_Read,0,0,width,height,BufferLandsat_0[k],width,height,GDT_Float32,0,0);
	}	

	GDALDataset *MODIS0 = (GDALDataset*) GDALOpen(BufferIn1,GA_ReadOnly);
	float** BufferModis_0 = new float*[BandNum];
	for( b=0;b<BandNum;b++)
	{
		BufferModis_0[b] = new float[width*height];	
	}
	
	for( k=0;k<BandNum;k++)
	{
		GDALRasterBand* hInBand1 = MODIS0->GetRasterBand(k+1);
		hInBand1->RasterIO(GF_Read,0,0,width,height,BufferModis_0[k],width,height,GDT_Float32,0,0);		
	}	

	GDALDataset *Landsat1 = (GDALDataset*) GDALOpen(BufferIn2,GA_ReadOnly);
	float** BufferLandsat_1 = new float*[BandNum];
	for( b=0;b<BandNum;b++)
	{
		BufferLandsat_1[b] = new float[width*height];	
	}
	
	for(k=0;k<BandNum;k++)
	{
		GDALRasterBand* hInBand1 = Landsat1->GetRasterBand(k+1);
		hInBand1->RasterIO(GF_Read,0,0,width,height,BufferLandsat_1[k],width,height,GDT_Float32,0,0);
	}	

	GDALDataset *MODIS1 = (GDALDataset*) GDALOpen(BufferIn3,GA_ReadOnly);
	float** BufferModis_1 = new float*[BandNum];
	for( b=0;b<BandNum;b++)
	{
		BufferModis_1[b] = new float[width*height];	
	}
	
	for( k=0;k<BandNum;k++)
	{
		GDALRasterBand* hInBand1 = MODIS1->GetRasterBand(k+1);
		hInBand1->RasterIO(GF_Read,0,0,width,height,BufferModis_1[k],width,height,GDT_Float32,0,0);		
	}	
	
	GDALDataset *MODIS2 = (GDALDataset*) GDALOpen(BufferIn4,GA_ReadOnly);
	
	float** BufferModis_2 = new float*[BandNum];
	for( b=0;b<BandNum;b++)
	{
		BufferModis_2[b] = new float[width*height];	
	}
	
	for( k=0;k<BandNum;k++)
	{
		GDALRasterBand* hInBand1 = MODIS2->GetRasterBand(k+1);
		hInBand1->RasterIO(GF_Read,0,0,width,height,BufferModis_2[k],width,height,GDT_Float32,0,0);		
		
	}
	
	GDALDataset *LandsatDs;
	char* driverName = "GTiff";
	GDALDriver *pDriver = (GDALDriver*)GDALGetDriverByName(driverName);
	LandsatDs = pDriver->Create(BufferOut,width,height,BandNum,GDT_Float64,NULL);
	double* geos=new double[6];
	Landsat0->GetGeoTransform(geos);
	LandsatDs->SetGeoTransform(geos);
	LandsatDs->SetProjection(Landsat0->GetProjectionRef());
	
	float** BufferOutColor = new float*[BandNum];
	for( b=0;b<BandNum;b++)
	{
		BufferOutColor[b] = new float[width*height];
	}
	//e.Blending2(BufferLandsat_0,BufferModis_0,BufferLandsat_1,BufferModis_1,BufferModis_2,BufferOutColor,height,width,win_size,flag, L_err, M_err, Para_h,BandNum,1.0);
	long now1 = clock();
	 runtest(BufferLandsat_0,BufferModis_0,BufferLandsat_1,BufferModis_1,BufferModis_2,BufferOutColor,height,width,win_size,M_err,BandNum,class_num,_nodata);
	  printf("GPU运行时间为：%dms\n", int(((double)(clock() - now1)) / CLOCKS_PER_SEC * 1000));
	for (b=0;b<BandNum;b++)
	{
		GDALRasterBand* HOut = LandsatDs->GetRasterBand(b+1);
		HOut->RasterIO(GF_Write,0,0,width,height,BufferOutColor[b],width,height,GDT_Float32,0,0);
	}
	GDALClose(Landsat0);
	GDALClose(MODIS0);
	GDALClose(Landsat1);
	GDALClose(MODIS1);
	GDALClose(MODIS2);
	GDALClose(LandsatDs);

	for (b=0;b<BandNum;b++)
	{
		delete []BufferLandsat_0[b];
		delete []BufferModis_0[b];
		delete []BufferLandsat_1[b];
		delete []BufferModis_1[b];
		delete []BufferModis_2[b];
		delete []BufferOutColor[b];
	}
	delete []BufferLandsat_0;
	delete [] BufferModis_0;
	delete []BufferLandsat_1;
	delete [] BufferModis_1;
	delete [] BufferModis_2;
	delete [] BufferOutColor;
}


