/*************************************************************************
 ����ʹ������ԭ�ӵı�׼����, ����ѧϰʹ��
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
void BEEP_GPIO2(void)//����������
{	
GPIO_InitTypeDef GPIO_InitStruct;
    // ʹ��GPIOB�˿�ʱ��
    RCC_APB2PeriphClockCmd(RCC_APB2Periph_GPIOB, ENABLE);

    // ����PB9����Ϊ���ģʽ������������Ϊ10MHz
    
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

uint32_t aun_ir_buffer[500]; 	//IR����,�����ԭʼ���ݻ�����
uint32_t aun_red_buffer[500];  	//Red����,���ԭʼ���ݻ�����
int32_t n_ir_buffer_length;    	//���ݳ���
	
int32_t n_sp02; 				//Ѫ��ֵ
int32_t n_heart_rate;   		//����ֵ
	
int8_t ch_spo2_valid;   		//�㷨ִ�гɹ���־λ
int8_t ch_hr_valid;    			

int main(void)
{ 
	uint32_t un_min, un_max; 	//ԭʼ���ݵ����/��Сֵ
	uint32_t Wave_Range;		//������ʾ���㱶��
	u16 i;						//
	
	s8 Wave_sum;				//OLED���δ�С
	u8 temp[6];					//ƴ�ֽ�
	u8 str[30];					//�ַ�����ʾ����
	u8 dis_hr=0,dis_spo2=0;
	u8 wave;
	
	NVIC_Configuration();     	//��������ʼ��
	delay_init();	    		  
	uart_init(115200);
	LED_Init();
	OLED_Init();
	BEEP_GPIO2();
	max30102_init();
	maxim_max30102_read_reg(REG_REV_ID,temp);//��ID,�ɲ�Ҫ
	
	sprintf((char *)str,"  HR:   SpO2:  ");
	OLED_ShowString(0,0,str,16);
	OLED_Refresh_Gram();		//������ʾ��OLED	 

	n_ir_buffer_length=500; 	//�������鳤��
	while(1)
	{
		i=0;					//���ó�ʼֵ
		un_max=0;
        un_min=PPG_MIN_COUNT;
        
        for(i=100;i<500;i++)	//������ǰ��100->Ŀ���Ƕ���500�������е�ǰ��100������,������100���²ɼ���ֵ
        {
            aun_red_buffer[i-100]=aun_red_buffer[i];
            aun_ir_buffer[i-100]=aun_ir_buffer[i];
            
			if(i > 400)
			{
				if(un_min>aun_ir_buffer[i])  
					un_min=aun_ir_buffer[i];	//�ҵ���Сֵ
				if(un_max<aun_ir_buffer[i])
					un_max=aun_ir_buffer[i];	//�ҵ����ֵ
			}
			
        }
		Wave_Range = un_max-un_min; //���㲨�岨�ȵĲ�ֵ
		Wave_Range = Wave_Range/14; //���㱶��

        for(i=400;i<500;i++)  	//���100�������ݵ�����
        {
			wave++;				//����OLEDˢ���ٶ�
            while(MAX30102_INT==1);
            max30102_FIFO_ReadBytes(REG_FIFO_DATA,temp);
			aun_ir_buffer[i] =  (long)((long)((long)temp[0]&0x03)<<16) | (long)temp[1]<<8 | (long)temp[2]; //�ɼ��ͺϲ��õ�ԭʼ���� 
			aun_red_buffer[i] = (long)((long)((long)temp[3] & 0x03)<<16) |(long)temp[4]<<8 | (long)temp[5];  

			if(aun_ir_buffer[i] > 5000) //��ָ�ſ�ʱһ�㶼����5000
			{
				if(aun_ir_buffer[i] > un_max) un_max = (aun_ir_buffer[i] + 100); //ˢ�±���
				else if(aun_ir_buffer[i] < un_min)
				{
					un_min = aun_ir_buffer[i] - 200;
					Wave_Range = un_max-un_min;
					Wave_Range = Wave_Range/14;
				}
				
				Wave_sum = (un_max - aun_ir_buffer[i])/Wave_Range; //�õ�OLED��ʾ�Ĳ���ֵ
			}
			else Wave_sum=1;
	
			if(Wave_sum>16) Wave_sum = 16; 							//��ֵ����,����
			printf("%d, %d, %d, %d \r\n",aun_ir_buffer[i],un_max,un_min,Wave_Range);
			if(wave==PPG_WAVE_COUNT)
			{
				OLED_wave(Wave_sum);  								//��ʾ����
				wave=0;
			}
		}
		//ִ���㷨,ȥֱ��,�˲�,���㲨�η�ֵ��
        maxim_heart_rate_and_oxygen_saturation(aun_ir_buffer, n_ir_buffer_length, aun_red_buffer, &n_sp02, &ch_spo2_valid, &n_heart_rate, &ch_hr_valid);
		
		if(ch_hr_valid==1 && ch_spo2_valid==1&& n_sp02<100&&n_sp02>70) //������ʾ����, ��ѧ������Ȼ��Ӱ��
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
		OLED_Refresh_Gram();//������ʾ��OLED	
	}
}
void OLED_wave(u8 Wave_sum)
{
		static u8 i;		  
	
		OLED_WR_Byte (0xb0,OLED_CMD);    						//����ҳ��ַ��0��
		OLED_WR_Byte ((i & 0x0f),OLED_CMD);      				//������ʾλ�á��е͵�ַ
		OLED_WR_Byte (((i & 0xf0) >> 4) |0x10,OLED_CMD);      	//������ʾλ�á��иߵ�ַ   
		OLED_WR_Byte(PPG_wave_fig[Wave_sum*2],OLED_DATA); 
	
		OLED_WR_Byte (0xb1,OLED_CMD);    						//����ҳ��ַ��1��
		OLED_WR_Byte ((i & 0x0f),OLED_CMD);      				//������ʾλ�á��е͵�ַ
		OLED_WR_Byte (((i & 0xf0) >> 4) |0x10,OLED_CMD);      	//������ʾλ�á��иߵ�ַ 
		OLED_WR_Byte(PPG_wave_fig[Wave_sum*2+1],OLED_DATA);		//��ʾ����,PPG_wave_fig[n]��ֵ��OLED�Ĵ������õ�ɨ�跽ʽ���
		i++;
		if(i>127) i=0;

}

