/*************************************************************************
 程序使用正点原子的标准库框架, 仅供学习使用
**************************************************************************************/
#include "led.h"
#include "delay.h"
#include "sys.h"
#include "usart.h"
#include "max30102.h" 
#include "myiic.h"
#include "algorithm.h"
#include "oled.h"
#include <stdlib.h>
//#include <stdio.h>

#define MAX_BRIGHTNESS 255
#define PPG_WAVE_COUNT 	4
#define PPG_MIN_COUNT 	0X3FFFF
#define BEEP_Set() GPIO_SetBits(GPIOB,GPIO_Pin_9)
#define BEEP_off() GPIO_ResetBits(GPIOB,GPIO_Pin_9)
void OLED_wave(u8 Wave_sum);
void BEEP_GPIO2(void)//蜂鸣器配置
{	
GPIO_InitTypeDef GPIO_InitStruct;
    // 使能GPIOB端口时钟
    RCC_APB2PeriphClockCmd(RCC_APB2Periph_GPIOB, ENABLE);

    // 配置PB9引脚为输出模式，最大输出速率为10MHz
    
    GPIO_InitStruct.GPIO_Mode = GPIO_Mode_Out_PP;
    GPIO_InitStruct.GPIO_Pin = GPIO_Pin_9;
    GPIO_InitStruct.GPIO_Speed = GPIO_Speed_10MHz;
    GPIO_Init(GPIOB, &GPIO_InitStruct);

   
}
const unsigned char PPG_wave_fig[34]= //0-16
{	0x00,0x00,
	0x01,0x00,
	0x03,0x00,
	0x07,0x00,
	0x0f,0x00,
	0x1f,0x00,
	0x3f,0x00,
	0x7f,0x00,
	0xff,0x00,
	0xff,0x01,
	0xff,0x03,
	0xff,0x07,
	0xff,0x0f,
	0xff,0x1f,
	0xff,0x3f,
	0xff,0x7f,
	0xff,0xff};

uint32_t aun_ir_buffer[500]; 	//IR数组,红外光原始数据缓存区
uint32_t aun_red_buffer[500];  	//Red数组,红光原始数据缓存区
int32_t n_ir_buffer_length;    	//数据长度
	
int32_t n_sp02; 				//血氧值
int32_t n_heart_rate;   		//心率值
	
int8_t ch_spo2_valid;   		//算法执行成功标志位
int8_t ch_hr_valid;    			

int main(void)
{ 
	uint32_t un_min, un_max; 	//原始数据的最大/最小值
	uint32_t Wave_Range;		//用于显示计算倍率
	u16 i;						//
	
	s8 Wave_sum;				//OLED波形大小
	u8 temp[6];					//拼字节
	u8 str[30];					//字符串显示变量
	u8 dis_hr=0,dis_spo2=0;
	u8 wave;
	
	NVIC_Configuration();     	//各函数初始化
	delay_init();	    		  
	uart_init(115200);
	LED_Init();
	OLED_Init();
	BEEP_GPIO2();
	max30102_init();
	maxim_max30102_read_reg(REG_REV_ID,temp);//读ID,可不要
	
	sprintf((char *)str,"  HR:   SpO2:  ");
	OLED_ShowString(0,0,str,16);
	OLED_Refresh_Gram();		//更新显示到OLED	 

	n_ir_buffer_length=500; 	//更新数组长度
	while(1)
	{
		i=0;					//重置初始值
		un_max=0;
        un_min=PPG_MIN_COUNT;
        
        for(i=100;i<500;i++)	//数组往前移100->目的是丢弃500个数组中的前面100个数组,并更新100个新采集的值
        {
            aun_red_buffer[i-100]=aun_red_buffer[i];
            aun_ir_buffer[i-100]=aun_ir_buffer[i];
            
			if(i > 400)
			{
				if(un_min>aun_ir_buffer[i])  
					un_min=aun_ir_buffer[i];	//找到最小值
				if(un_max<aun_ir_buffer[i])
					un_max=aun_ir_buffer[i];	//找到最大值
			}
			
        }
		Wave_Range = un_max-un_min; //计算波峰波谷的差值
		Wave_Range = Wave_Range/14; //计算倍率

        for(i=400;i<500;i++)  	//填充100个新数据到数组
        {
			wave++;				//减慢OLED刷新速度
            while(MAX30102_INT==1);
            max30102_FIFO_ReadBytes(REG_FIFO_DATA,temp);
			aun_ir_buffer[i] =  (long)((long)((long)temp[0]&0x03)<<16) | (long)temp[1]<<8 | (long)temp[2]; //采集和合并得到原始数据 
			aun_red_buffer[i] = (long)((long)((long)temp[3] & 0x03)<<16) |(long)temp[4]<<8 | (long)temp[5];  

			if(aun_ir_buffer[i] > 5000) //手指放开时一般都大于5000
			{
				if(aun_ir_buffer[i] > un_max) un_max = (aun_ir_buffer[i] + 100); //刷新倍率
				else if(aun_ir_buffer[i] < un_min)
				{
					un_min = aun_ir_buffer[i] - 200;
					Wave_Range = un_max-un_min;
					Wave_Range = Wave_Range/14;
				}
				
				Wave_sum = (un_max - aun_ir_buffer[i])/Wave_Range; //得到OLED显示的波形值
			}
			else Wave_sum=1;
	
			if(Wave_sum>16) Wave_sum = 16; 							//数值过大,纠正
			printf("%d, %d, %d, %d \r\n",aun_ir_buffer[i],un_max,un_min,Wave_Range);
			if(wave==PPG_WAVE_COUNT)
			{
				OLED_wave(Wave_sum);  								//显示波形
				wave=0;
			}
		}
		//执行算法,去直流,滤波,计算波形幅值等
        maxim_heart_rate_and_oxygen_saturation(aun_ir_buffer, n_ir_buffer_length, aun_red_buffer, &n_sp02, &ch_spo2_valid, &n_heart_rate, &ch_hr_valid);
		
		if(ch_hr_valid==1 && ch_spo2_valid==1&& n_sp02<100&&n_sp02>70) //划定显示区间, 光学会受自然光影响
		{
			dis_hr = n_heart_rate;
			dis_spo2 = n_sp02;
		}
		else
		{
			dis_hr = 0;
			dis_spo2 = 0;
		}
		if(dis_hr == 0 || dis_spo2 == 0)  
		{
			sprintf((char *)str," --   -- ");
			OLED_ShowString(0,20,str,24);
		}
		else{
			sprintf((char *)str,"%3d  %3d",dis_hr,dis_spo2);
			OLED_ShowString(0,20,str,24);
		}
		if(dis_hr<60||dis_spo2>100||dis_spo2<94)
		BEEP_Set();
	else if(dis_hr>=60&&dis_hr<=100&&dis_spo2>=94)
		BEEP_off();
		OLED_Refresh_Gram();//更新显示到OLED	
	}
}
void OLED_wave(u8 Wave_sum)
{
		static u8 i;		  
	
		OLED_WR_Byte (0xb0,OLED_CMD);    						//设置页地址（0）
		OLED_WR_Byte ((i & 0x0f),OLED_CMD);      				//设置显示位置—列低地址
		OLED_WR_Byte (((i & 0xf0) >> 4) |0x10,OLED_CMD);      	//设置显示位置—列高地址   
		OLED_WR_Byte(PPG_wave_fig[Wave_sum*2],OLED_DATA); 
	
		OLED_WR_Byte (0xb1,OLED_CMD);    						//设置页地址（1）
		OLED_WR_Byte ((i & 0x0f),OLED_CMD);      				//设置显示位置—列低地址
		OLED_WR_Byte (((i & 0xf0) >> 4) |0x10,OLED_CMD);      	//设置显示位置—列高地址 
		OLED_WR_Byte(PPG_wave_fig[Wave_sum*2+1],OLED_DATA);		//显示数据,PPG_wave_fig[n]的值跟OLED寄存器配置的扫描方式相关
		i++;
		if(i>127) i=0;

}

