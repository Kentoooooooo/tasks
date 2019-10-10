#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <windows.h>
#include "dSFMT.h"

#define C 2.99792458e8
#define h 1e-10
#define dh 1e-4
#define EPS 1e-8
#define WGSD 1.0e4
#define OUT_NODE_STATE 0

struct parameter {
	double tau_L, tau_H, delay[2], beta, kappa, sigma, phi, epsilon;
};

struct LaserEOfeedback {
	struct parameter p;
	double x, y, *xdelay;
	int descretized_num[2], delay_index[2], long_delay;
};

struct Reservoir {
	int trainingNum, testNum, mask_num, node_num, allNum, mask_interval, node_interval, input_interval, mask_idx, node_idx, input_idx;
	double gamma, *mask, **train_node, **test_node, *target_train, *target_test, *weight, *input_train, *input_test, *output_train, *output_test, input_signal, input_MZM, input_MZM_phi, NMSE_train, NMSE_test;
};

struct Histdata {
	double min, max, div, minData, maxData, median, average, variance, entropy;
	int bin, *data, n;
};

struct cc_data {
	double sum[2], sq_sum[2], product_sum, correlation, difference;
	int num;
};

void CalcLaserEOfeedback(struct LaserEOfeedback *laser, double input_signal);
void SetRC_Param(struct Reservoir *rc, double mask_interval, int mask_num, int mask_node_ratio);
void MakeMask(struct Reservoir *rc);
void ReadTWFM(struct Reservoir *rc);
void QR(struct Reservoir *rc);
void EvaluateReservoirTrain(struct Reservoir *rc, char *filename);
void EvaluateReservoirTest(struct Reservoir *rc, char *filename);
void Surrogate_Node(struct Reservoir *rc);
void free_rc_memory(struct Reservoir *rc);
void SetHistgramParameter(struct Histdata *histgram, double min, double max, int bin);
void InputHistgram(struct Histdata *histgram, double inputData);
void OutputHistgram(struct Histdata *histgram, char *filename);
void Set_CC_fixed_shift(struct cc_data *cc);
void Input_CC_fixed_shift(struct cc_data *cc, double input1, double input2);
void Output_CC_fixed_shift(struct cc_data *cc);

int main(int argc, char *argv[]) {
	struct LaserEOfeedback laser[2];
	struct Reservoir rc;
	struct Histdata hist[4];
	struct cc_data cc[2];
	double param, param_min, param_max, param_interval, t_transient, t_plot;
	double buff[8], *wfmData[8], node_vc_difference;
	int i, j, k, n, p, s, t, transient, calc_num, divi, out, para_flag, iteration, plotnum;
	char buffer_ch[128], filename[128];
	clock_t start, end;
	FILE *inputfp, *outputfp[4];
	
	//	ファイルからのデータの読み込み
	if (argc == 2) {
		inputfp = fopen(argv[1], "r");
		if(inputfp == NULL) {
			fprintf(stderr, "%sを開くことができません", argv[1]);
			exit(2);
		}
		fscanf(inputfp, "%lf", &(laser[0].p.tau_L));
		fscanf(inputfp, "%lf", &(laser[0].p.tau_H));
		fscanf(inputfp, "%lf %lf %lf %lf", &(laser[0].p.delay[0]), &(laser[1].p.delay[0]), &(laser[0].p.delay[1]), &(laser[1].p.delay[1]));
		fscanf(inputfp, "%lf", &(laser[0].p.beta));
		fscanf(inputfp, "%lf %lf", &(laser[0].p.kappa), &(laser[1].p.kappa));
		fscanf(inputfp, "%lf %lf", &(laser[0].p.sigma), &(laser[1].p.sigma));
		fscanf(inputfp, "%lf %lf", &(laser[0].p.phi), &(laser[1].p.phi));
		fscanf(inputfp, "%lf %lf %lf", &param_min, &param_max, &param_interval);
		fscanf(inputfp, "%lf %lf", &t_transient, &t_plot);
		fscanf(inputfp, "%d", &divi);
		
		fclose(inputfp);
	} else {
		fprintf(stderr, "引数のエラーです\n");
		exit(1);
	}
	
	const int amp_max_int = (param_max * param_interval >= 0.0) ? (int)(param_max / param_interval + EPS) : (int)(param_max / param_interval - EPS);
	const int amp_min_int = (param_min * param_interval >= 0.0) ? (int)(param_min / param_interval + EPS) : (int)(param_min / param_interval - EPS);
	
	//	出力モード選択部分
	printf("Select a calculation-mode (1 or 2)\n");
	printf("\t 1 --> Temporal waveforms\n");
	printf("\t 2 --> Parameter change\n");
	printf("inputfp --> ");
	scanf("%d", &out);
	if(out == 1) {
		printf("Select an output-mode (1 or 2 or 3)\n");
		printf("\t 1 --> Temporal waveforms of x(t) and y(t)\n");
		printf("\t 2 --> Reservoir Computing\n");
		printf("inputfp --> ");
		scanf("%d", &para_flag);
		if(para_flag == 1) {
			outputfp[0] = fopen("twfm.txt", "w");
			fprintf(outputfp[0], "Time [\\fm\\ns]\tx\\d1\\n(t)\tx\\d2\\n(t)\tMask\tIdeal\tInput\tInput MZM output\tNode x\\d1\\n(t)\tNode x\\d2\\n(t)\n");
		} else if(para_flag == 2) {
		} else {
			fprintf(stderr, "error : para_flag = %d\n", para_flag);
			exit(1);
		}
	} else if(out == 2) {
		printf("Select a parameter which you change (1 or 2 or 3)\n");
		printf("\t 1 --> Feedback strength beta\n");
		printf("\t 2 --> Feedback and coupling delay time (all)\n");
		printf("\t 3 --> Gamma\n");
		printf("\t 4 --> Kappa (all)\n");
		printf("\t 5 --> Sigma for both MZMs\n");
		printf("\t 6 --> Phi_0 for MZM2\n");
		printf("\t 7 --> Num of nodes\n");
		printf("\t 8 --> Node interval\n");
		printf("\t 9 --> Feedback and coupling delay time (MZM1)\n");
		printf("\t 10 --> Feedback and coupling delay time (MZM2)\n");
		printf("\t 11 --> Feedback delay time (MZM2)\n");
		printf("\t 12 --> Mask pattern\n");
		printf("inputfp --> ");
		scanf("%d", &para_flag);
		if(para_flag == 1) sprintf(buffer_ch, "Feedback strength \\fb\\n");
		if(para_flag == 2) sprintf(buffer_ch, "Feedback and coupling delay time (all) \\ft\\n\\d1\\n [ms]");
		if(para_flag == 3) sprintf(buffer_ch, "\\fg\\n");
		if(para_flag == 4) sprintf(buffer_ch, "\\fk\\n for all");
		if(para_flag == 5) sprintf(buffer_ch, "\\fs\\n for both MZMs");
		if(para_flag == 6) sprintf(buffer_ch, "\\fp\\n\\d0\\n for MZM2");
		if(para_flag == 7) sprintf(buffer_ch, "Number of nodes");
		if(para_flag == 8) sprintf(buffer_ch, "Node interval [\\fm\\ns]");
		if(para_flag == 9) sprintf(buffer_ch, "Delay time (MZM1) \\ft\\n\\d1\\n [ms]");
		if(para_flag == 10) sprintf(buffer_ch, "Delay time (MZM2) \\ft\\n\\d2\\n [ms]");
		if(para_flag == 11) sprintf(buffer_ch, "Delay time (MZM2) \\ft\\n\\d2\\n [ms]");
		if(para_flag == 12) sprintf(buffer_ch, "Mask pattern");
		outputfp[0] = fopen("param_change.txt", "w");
		fprintf(outputfp[0], "%s\tNMSE for training\tNMSE for test\tVariance of MZM1 nodes\t\tVariance of MZM2 nodes\t\tVariance of all nodes\tcc between twfms\tcc between nodes\tDifference between twfms\tDifference between nodes\tDifference between node vectors\n", buffer_ch);
		fclose(outputfp[0]);
		
		if(OUT_NODE_STATE == 1) CreateDirectory(".\\node_state", NULL);
	} else {
		fprintf(stderr, "error : out = %d\n", out);
		exit(1);
	}
	
	/* ここから初期化部分 */
	laser[0].p.tau_L = 1.0 / (2.0 * M_PI * laser[0].p.tau_L);
	laser[0].p.tau_H = 1.0 / (2.0 * M_PI * laser[0].p.tau_H);
	laser[1].p.tau_L = laser[0].p.tau_L;
	laser[1].p.tau_H = laser[0].p.tau_H;
	laser[1].p.beta = laser[0].p.beta;
	laser[0].p.phi = laser[0].p.phi * M_PI;
	laser[1].p.phi = laser[1].p.phi * M_PI;
	laser[0].p.epsilon = laser[1].p.epsilon = laser[0].p.tau_L / laser[0].p.tau_H;
	printf("epsilon = %e, beta = %e\n", laser[0].p.epsilon, laser[0].p.beta);
	
	for(i = 0; i < 2; i++) {
		if(laser[i].p.delay[0] >= laser[i].p.delay[1]) laser[i].long_delay = 0;
		else laser[i].long_delay = 1;
		for(j = 0; j < 2; j++) laser[i].descretized_num[j] = (int)(laser[i].p.delay[j] / h + EPS);
		laser[i].xdelay = (double*)malloc(sizeof(double) * laser[i].descretized_num[laser[i].long_delay]);
		printf("Laser %d\nDelay time = %e (Descretized num = %d)\n", i, laser[i].p.delay[0], laser[i].descretized_num[0]);
	}
	
	SetRC_Param(&rc, 0.1e-6, 400, 1);//mask interval, number of mask, ratio of node to mask
	
	transient = (int)(t_transient / dh + EPS);
	plotnum = (int)(t_plot / dh + EPS);
	printf("Number of node: %d\n", rc.node_num);
	printf("Number of one mask point: %d\n", rc.node_interval);
	/* ここまで初期化部分 */
	
	//	メイン計算
	for(s = amp_min_int; s <= amp_max_int; s++) {
		param = param_interval * (double)s;
		
		start = clock();
		
		if(out == 2) {
			if(para_flag == 1) {
				laser[0].p.beta = laser[1].p.beta = param;
				printf("\nbeta = %e\n", laser[0].p.beta);
			}
			if(para_flag == 2) {
				laser[0].p.delay[0] = laser[0].p.delay[1] = laser[1].p.delay[0] = laser[1].p.delay[1] = param;
				for(i = 0; i < 2; i++) {
					laser[i].long_delay = 0;
					for(j = 0; j < 2; j++) laser[i].descretized_num[j] = (int)(laser[i].p.delay[j] / h + EPS);
					free(laser[i].xdelay);
					laser[i].xdelay = (double*)malloc(sizeof(double) * laser[i].descretized_num[laser[i].long_delay]);
				}
				printf("\nDelay time of all = %e (Descretized num = %d)\n", laser[0].p.delay[0], laser[0].descretized_num[0]);
				param *= dh / h;
			}
			if(para_flag == 3) {
				rc.gamma = param / 255.0;
				printf("\ngamma = %e\n", param);
			}
			if(para_flag == 4) {
				laser[0].p.kappa = param;
				// laser[1].p.kappa = param;
				printf("\nkappa_2 (kappa_1 = %e) = %e\n", laser[1].p.kappa, laser[0].p.kappa);
			}
			if(para_flag == 5) {
				laser[1].p.sigma = param;
				laser[0].p.sigma = param;
				printf("\nsigma = %e\n", laser[0].p.sigma);
			}
			if(para_flag == 6) {
				laser[1].p.phi = M_PI * param;
				printf("\nphi_0 for MZM2 = %e * PI\n", param);
			}
			if(para_flag == 7) {
				free_rc_memory(&rc);
				// SetRC_Param(&rc, 10.0e-6 / param, (int)(param));
				SetRC_Param(&rc, 0.4e-6, param, 2);
				printf("\nNum of nodes = %d (node interval = %e)\n", rc.node_num, rc.node_interval * h);
			}
			if(para_flag == 8) {
				free_rc_memory(&rc);
				SetRC_Param(&rc, 0.4e-6, 25, 2);
				printf("\nNode interval = %e\n", rc.node_interval * h);
				param *= 1.0e6;
			}
			if(para_flag == 9) {
				laser[0].p.delay[0] = laser[0].p.delay[1] = param;
				laser[0].long_delay = 0;
				for(j = 0; j < 2; j++) laser[0].descretized_num[j] = (int)(laser[0].p.delay[j] / h + EPS);
				free(laser[0].xdelay);
				laser[0].xdelay = (double*)malloc(sizeof(double) * laser[0].descretized_num[laser[0].long_delay]);
				printf("\nDelay time for MZM1 = %e (Descretized num = %d)\n", laser[0].p.delay[0], laser[0].descretized_num[0]);
				param *= dh / h;
			}
			if(para_flag == 10) {
				laser[1].p.delay[0] = laser[1].p.delay[1] = param;
				laser[1].long_delay = 0;
				for(j = 0; j < 2; j++) laser[1].descretized_num[j] = (int)(laser[1].p.delay[j] / h + EPS);
				free(laser[1].xdelay);
				laser[1].xdelay = (double*)malloc(sizeof(double) * laser[1].descretized_num[laser[1].long_delay]);
				printf("\nDelay time for MZM2 = %e (Descretized num = %d)\n", laser[1].p.delay[0], laser[1].descretized_num[0]);
				param *= dh / h;
			}
			if(para_flag == 11) {
				laser[0].p.delay[1] = laser[1].p.delay[0] = param;
				for(i = 0; i < 2; i++) {
					if(laser[i].p.delay[0] >= laser[i].p.delay[1]) laser[i].long_delay = 0;
					else laser[i].long_delay = 1;
					for(j = 0; j < 2; j++) laser[i].descretized_num[j] = (int)(laser[i].p.delay[j] / h + EPS);
					free(laser[i].xdelay);
					laser[i].xdelay = (double*)malloc(sizeof(double) * laser[i].descretized_num[laser[i].long_delay]);
				}
				printf("\nPropagation time from MZM1 = %e (Descretized num = %d)\n", laser[0].p.delay[1], laser[0].descretized_num[1]);
				printf("Delay time for MZM2 = %e (Descretized num = %d)\n", laser[1].p.delay[0], laser[1].descretized_num[0]);
				param *= dh / h;
			}
			if(para_flag == 12) {
				printf("\nMask pattern = %e\n", param);
				for(i = 0, buff[0] = param; i < rc.mask_num; i++) {
					if((int)buff[0] % 2 == 1) rc.mask[i] = 0.1;
					else rc.mask[i] = -0.1;
					buff[0] = (int)(buff[0] / 2);
					printf("%+1.1f, ", rc.mask[i]);
				}
				printf("\n");
			}
		}
		
		for(i = 0; i < 2; i++) {
			laser[i].x = 0.05 + (double)i * 0.05;
			laser[i].y = 0.0;
			for(j = 0; j < laser[i].descretized_num[laser[i].long_delay]; j++) laser[i].xdelay[j] = 0.05 + (double)i * 0.05;
			for(j = 0; j < 2; j++) laser[i].delay_index[j] = laser[i].descretized_num[laser[i].long_delay] - laser[i].descretized_num[j];
		}
		init_gen_rand(222);
		
		SetHistgramParameter(&hist[0], -1.0, 2.0, 100);
		SetHistgramParameter(&hist[1], -1.0, 2.0, 100);
		SetHistgramParameter(&hist[2], -1.0, 2.0, 100);
		SetHistgramParameter(&hist[3], 0.0, 5.0, 50);
		
		Set_CC_fixed_shift(&cc[0]);
		Set_CC_fixed_shift(&cc[1]);
		
		iteration = -10;
		for(i = 0; i < 8; i++) buff[i] = 0.0;
		
		for(i = 0; i < transient; i++) CalcLaserEOfeedback(laser, 0.5);
		
		calc_num = rc.trainingNum * rc.input_interval;
		rc.mask_idx = rc.node_idx = rc.input_idx = 0;
		
		for(i = 1; i <= calc_num; i++) {
			
			rc.input_signal = rc.gamma * rc.mask[rc.mask_idx] * rc.input_train[rc.input_idx];
			rc.input_MZM = cos(rc.input_signal * M_PI / 4.0 + rc.input_MZM_phi);
			rc.input_MZM *= rc.input_MZM;
			
			// CalcLaserEOfeedback(laser, 0.5);
			CalcLaserEOfeedback(laser, rc.input_MZM);
			
			Input_CC_fixed_shift(&cc[0], laser[0].x, laser[1].x);
			
			if(out == 1) {
				if(para_flag == 1) {
					if(i % divi == 0) {
						// if(laser[0].p.delay[0] <= laser[1].p.delay[0]) buff[0] = laser[0].x + laser[1].xdelay[(laser[1].delay_index[0] + laser[0].descretized_num[0]) % laser[1].descretized_num[laser[1].long_delay]];
						// else buff[0] = laser[0].xdelay[(laser[0].delay_index[0] + laser[1].descretized_num[0]) % laser[0].descretized_num[laser[0].long_delay]] + laser[1].x;
						if(i % rc.node_interval == rc.node_interval / 2) 
							fprintf(outputfp[0], "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", i * dh, laser[0].x, laser[1].x, rc.mask[rc.mask_idx], rc.input_train[rc.input_idx], rc.input_signal, rc.input_MZM, laser[0].x, laser[1].x);
						else fprintf(outputfp[0], "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", i * dh, laser[0].x, laser[1].x, rc.mask[rc.mask_idx], rc.input_train[rc.input_idx], rc.input_signal, rc.input_MZM);
					}
					if(i == plotnum) {
						fclose(outputfp[0]);
						free_rc_memory(&rc);
						for(i = 0; i < 2; i++) free(laser[i].xdelay);
						return 0;
					}
				}
			}
			
			if(i % rc.node_interval == rc.node_interval / 2) {
				rc.train_node[rc.input_idx][rc.node_idx] = laser[0].x;
				rc.train_node[rc.input_idx][rc.node_num + rc.node_idx] = laser[1].x;
				rc.node_idx = (rc.node_idx + 1) % rc.node_num;
				InputHistgram(&hist[0], laser[0].x);
				InputHistgram(&hist[1], laser[1].x);
				InputHistgram(&hist[2], laser[0].x);
				InputHistgram(&hist[2], laser[1].x);
				Input_CC_fixed_shift(&cc[1], laser[0].x, laser[1].x);
			}
			
			if(i % rc.mask_interval == 0) rc.mask_idx = (rc.mask_idx + 1) % rc.mask_num;
			if(i % rc.input_interval == 0) {
				for(j = 0, node_vc_difference = 0.0; j < rc.node_num; j++) node_vc_difference += pow(rc.train_node[rc.input_idx][j] - rc.train_node[rc.input_idx][rc.node_num + j], 2);
				InputHistgram(&hist[3], sqrt(node_vc_difference));
				rc.input_idx++;
			}
		}
		
		QR(&rc);
		if(out == 1) EvaluateReservoirTrain(&rc, "rc_training_result.txt\0");
		else EvaluateReservoirTrain(&rc, "0");
		
		calc_num = rc.testNum * rc.input_interval;
		rc.mask_idx = rc.node_idx = rc.input_idx = 0;
		
		for(i = 1; i <= calc_num; i++) {
			
			rc.input_signal = rc.gamma * rc.mask[rc.mask_idx] * rc.input_test[rc.input_idx];
			rc.input_MZM = cos(rc.input_signal * M_PI / 4.0 + rc.input_MZM_phi);
			rc.input_MZM *= rc.input_MZM;
			
			// CalcLaserEOfeedback(laser, 0.5);
			CalcLaserEOfeedback(laser, rc.input_MZM);
			
			Input_CC_fixed_shift(&cc[0], laser[0].x, laser[1].x);
			
			if(i % rc.node_interval == rc.node_interval / 2) {
				rc.test_node[rc.input_idx][rc.node_idx] = laser[0].x;
				rc.test_node[rc.input_idx][rc.node_num + rc.node_idx] = laser[1].x;
				rc.node_idx = (rc.node_idx + 1) % rc.node_num;
				InputHistgram(&hist[0], laser[0].x);
				InputHistgram(&hist[1], laser[1].x);
				InputHistgram(&hist[2], laser[0].x);
				InputHistgram(&hist[2], laser[1].x);
				Input_CC_fixed_shift(&cc[1], laser[0].x, laser[1].x);
			}
			
			if(i % rc.mask_interval == 0) rc.mask_idx = (rc.mask_idx + 1) % rc.mask_num;
			if(i % rc.input_interval == 0) {
				for(j = 0, node_vc_difference = 0.0; j < rc.node_num; j++) node_vc_difference += pow(rc.test_node[rc.input_idx][j] - rc.test_node[rc.input_idx][rc.node_num + j], 2);
				InputHistgram(&hist[3], sqrt(node_vc_difference));
				for(j = 0, rc.output_test[rc.input_idx] = 0.0; j < rc.node_num; j++) rc.output_test[rc.input_idx] += rc.weight[j] * rc.test_node[rc.input_idx][j] + rc.weight[j + rc.node_num] * rc.test_node[rc.input_idx][j + rc.node_num];
				rc.input_test[rc.input_idx + 1] = rc.output_test[rc.input_idx];
				rc.input_idx++;
			}
		}
		
		Output_CC_fixed_shift(&cc[0]);
		Output_CC_fixed_shift(&cc[1]);
		printf("===  CC  === \t%e\t%e\n", cc[0].correlation, cc[1].correlation);
		printf("=== DIFF === \t%e\t%e\n", cc[0].difference, cc[1].difference);
		
		end = clock();
		printf("%f\n", (double)(end - start) / CLOCKS_PER_SEC);
		
		// Surrogate_Node(&rc);
		
		if(out == 1) {
			if(para_flag == 2) {
				// outputfp[0] = fopen("node_state.txt", "w");
				// fprintf(outputfp[0], "Index\tNode\n");
				// for(i = 0; i < rc.allNum; i++) for(j = 0; j < rc.node_num; j++) fprintf(outputfp[0], "%d\t%e\n", i * rc.node_num + j, rc.nodeState[i][j]);
				// fclose(outputfp[0]);
				
				OutputHistgram(&hist[0], "histgram_laser1_node.txt\0");
				OutputHistgram(&hist[1], "histgram_laser2_node.txt\0");
				OutputHistgram(&hist[2], "histgram_all_node.txt\0");
				OutputHistgram(&hist[3], "histgram_node_vc_difference.txt\0");
				
				EvaluateReservoirTest(&rc, "rc_test_result.txt\0");
			}
			break;
		}
		if(out == 2) {
			OutputHistgram(&hist[0], "0");
			OutputHistgram(&hist[1], "0");
			OutputHistgram(&hist[2], "0");
			OutputHistgram(&hist[3], "0");
			
			EvaluateReservoirTest(&rc, "rc_test_result.txt\0");
			outputfp[0] = fopen("param_change.txt", "r+");
			fseek(outputfp[0], 0, SEEK_END);
			fprintf(outputfp[0], "%e\t%e\t%e\t", param, rc.NMSE_train, rc.NMSE_test);
			fprintf(outputfp[0], "%e\t%e\t%e\t", hist[0].variance, hist[1].variance, hist[2].variance);
			fprintf(outputfp[0], "%e\t%e\t%e\t%e\t%e\n", cc[0].correlation, cc[1].correlation, cc[0].difference, cc[1].difference, hist[3].average);
			fclose(outputfp[0]);
			
			if(OUT_NODE_STATE == 1) {
				sprintf(filename, ".\\node_state\\node_state_%03d.txt\0", s - amp_min_int);
				outputfp[1] = fopen(filename, "w");
				for(i = 0; i < rc.trainingNum; i++) {
					for(j = 0; j < rc.node_num * 2; j++) fprintf(outputfp[1], "%e\t", rc.train_node[i][j]);
					fprintf(outputfp[1], "\n");
				}
				fclose(outputfp[1]);
			}
		}
	}
	
	free_rc_memory(&rc);
	
	for(i = 0; i < 2; i++) free(laser[i].xdelay);
	
	return 0;
}

void free_rc_memory(struct Reservoir *rc) {
	int i;
	
	free(rc->target_train);
	free(rc->target_test);
	free(rc->input_train);
	free(rc->input_test);
	free(rc->output_train);
	free(rc->output_test);
	free(rc->mask);
	free(rc->weight);
	for(i = 0; i < rc->trainingNum; i++) free(rc->train_node[i]);
	free(rc->train_node);
	for(i = 0; i < rc->testNum; i++) free(rc->test_node[i]);
	free(rc->test_node);
	
	return;
}

void SetRC_Param(struct Reservoir *rc, double mask_interval, int mask_num, int mask_node_ratio) {
	int i;
	
	rc->allNum = 4000;
	rc->trainingNum = 3000;
	rc->testNum = rc->allNum - rc->trainingNum;
	rc->gamma = 1.0;//1.0 / 255.0;
	rc->input_MZM_phi = - M_PI / 4.0;
	
	// SantaFe
	rc->target_train = (double*)malloc(sizeof(double) * rc->trainingNum);
	rc->target_test = (double*)malloc(sizeof(double) * rc->testNum);
	rc->input_train = (double*)malloc(sizeof(double) * rc->trainingNum);
	rc->input_test = (double*)malloc(sizeof(double) * rc->testNum);
	rc->output_train = (double*)malloc(sizeof(double) * rc->trainingNum);
	rc->output_test = (double*)malloc(sizeof(double) * rc->testNum);
	ReadTWFM(rc);
	
	rc->mask_num = mask_num;
	rc->mask_interval = (int)(mask_interval / h + EPS);
	rc->input_interval = rc->mask_interval * rc->mask_num;
	
	rc->node_num = rc->mask_num * mask_node_ratio;
	rc->node_interval = rc->mask_interval / mask_node_ratio;
	
	// Make a mask for reservoir computing
	rc->mask = (double*)malloc(sizeof(double) * rc->mask_num);
	MakeMask(rc);
	
	rc->train_node = (double**)malloc(sizeof(double*) * rc->trainingNum);
	for(i = 0; i <  rc->trainingNum; i++) rc->train_node[i] = (double*)malloc(sizeof(double) * rc->node_num * 2);
	rc->test_node = (double**)malloc(sizeof(double*) * rc->testNum);
	for(i = 0; i <  rc->testNum; i++) rc->test_node[i] = (double*)malloc(sizeof(double) * rc->node_num * 2);
	rc->weight = (double*)malloc(sizeof(double) * rc->node_num * 2);
	
	return;
}

void MakeMask(struct Reservoir *rc) {
	int i, rand_num;
	
	init_gen_rand(1234);
	for(i = 0; i < rc->mask_num; i++) {
		// rand_num = (int)(11.0 * genrand_open_open());
		// if(rand_num == 0) rc->mask[i] = -1.0;
		// if(rand_num == 1) rc->mask[i] = -0.8;
		// if(rand_num == 2) rc->mask[i] = -0.6;
		// if(rand_num == 3) rc->mask[i] = -0.4;
		// if(rand_num == 4) rc->mask[i] = -0.2;
		// if(rand_num == 5) rc->mask[i] = -0.0;
		// if(rand_num == 6) rc->mask[i] = 0.2;
		// if(rand_num == 7) rc->mask[i] = 0.4;
		// if(rand_num == 8) rc->mask[i] = 0.6;
		// if(rand_num == 9) rc->mask[i] = 0.8;
		// if(rand_num == 10) rc->mask[i] = 1.0;
		
		// rand_num = (int)(6.0 * genrand_open_open());
		// if(rand_num == 0) rc->mask[i] = -1.0;
		// if(rand_num == 1) rc->mask[i] = -0.6;
		// if(rand_num == 2) rc->mask[i] = -0.2;
		// if(rand_num == 3) rc->mask[i] = -0.2;
		// if(rand_num == 4) rc->mask[i] = 0.6;
		// if(rand_num == 5) rc->mask[i] = 1.0;
		
		rand_num = (int)(4.0 * genrand_open_open());
		if(rand_num == 0) rc->mask[i] = -1.0;
		if(rand_num == 1) rc->mask[i] = -0.3;
		if(rand_num == 2) rc->mask[i] = 0.3;
		if(rand_num == 3) rc->mask[i] = 1.0;
		
		// rc->mask[i] = genrand_open_open();
		// if(rc->mask[i] >= 0.5){
			// rc->mask[i] = 1.0;
		// } else rc->mask[i] = -1.0;
	}
	
	return;
}

void ReadTWFM(struct Reservoir *rc) {
	int i, count, ret;
	double data, tmp[5000];
	FILE *inputfp;
	
	// inputfp = fopen("santafe.txt", "r");
	// inputfp = fopen("mg_twfm.txt", "r");
	inputfp = fopen("Lorenz.txt", "r");
	if(inputfp == NULL) {
		fprintf(stderr, "Can not open the file of santafe.txt.");
		exit(2);
	}
	
	count = 0;
	while((ret = fscanf(inputfp, "%lf", &data)) != EOF) {
		tmp[count] = (double)(data + 20) / 40;
		count++;
	}
	
	rc->target_train[2999] = tmp[3000];
	rc->input_train[0] = tmp[0];
	for(i = 1; i < 3000; i++) {
		rc->target_train[i - 1] = tmp[i];
		rc->input_train[i] = tmp[i];
	}
	
	rc->target_test[999] = tmp[4000];
	rc->input_test[0] = tmp[3000];
	for(i = 1; i < 1000; i++) {
		rc->target_test[i - 1] = tmp[i + 3000];
		rc->input_test[i] = tmp[i + 3000];
	}
	
	fclose(inputfp);
	
	return;
}

void QR(struct Reservoir *rc) {
	int i, j, k, sum_node = rc->node_num * 2;
	double sum, *yx, ridge = 0.001, *R[sum_node], *Q[sum_node], *QT[sum_node];
	
	for(i = 0; i < sum_node; i++) R[i] = (double*)malloc(sizeof(double) * sum_node);
	for(i = 0; i < sum_node; i++) Q[i] = (double*)malloc(sizeof(double) * sum_node);
	for(i = 0; i < sum_node; i++) QT[i] = (double*)malloc(sizeof(double) * sum_node);
	yx = (double*)malloc(sizeof(double) * sum_node);
	
	for(i = 0; i < sum_node; i++) {
		for(j = 0; j < sum_node; j++) {
			Q[i][j] = 0.0;
			for(k = 0; k < rc->trainingNum; k++) Q[i][j] += rc->train_node[k][i] * rc->train_node[k][j];
			Q[i][j] /= (double)rc->trainingNum;
			R[i][j] = 0.0;
		}
		yx[i] = 0.0;
		for(j = 0; j < rc->trainingNum; j++) yx[i] += rc->target_train[j] * rc->train_node[j][i];
		yx[i] /= (double)rc->trainingNum;
	}
	
	// for(i = 0; i < sum_node; i++) Q[i][i] += ridge;/* For ridge regression */
	
	for(i = 0; i < sum_node; i++) {
		for(j = 0; j < i; j++) {
			for(k = 0, R[j][i] = 0.0; k < sum_node; k++) R[j][i] += Q[k][i] * Q[k][j];
			for(k = 0; k < sum_node; k++) Q[k][i] -= Q[k][j] * R[j][i];
		}
		for(j = 0, R[i][i] = 0.0; j < sum_node; j++) R[i][i] += Q[j][i] * Q[j][i];
		R[i][i] = sqrt(R[i][i]);
		for(j = 0; j < sum_node; j++) Q[j][i] /= R[i][i];
	}
	
	for(i = 0; i < sum_node; i++) for(j = 0; j < sum_node; j++) QT[i][j] = Q[j][i];
	for(i = sum_node - 1; i >= 0; i--) {
		for(j = 0, sum = 0.0; j < sum_node; j++) sum += QT[i][j] * yx[j];
		for(j = i + 1; j < sum_node; j++) sum -= R[i][j] * rc->weight[j];
		rc->weight[i] = sum / R[i][i];
	}
	
	for(i = 0; i < sum_node; i++) free(R[i]);
	for(i = 0; i < sum_node; i++) free(Q[i]);
	for(i = 0; i < sum_node; i++) free(QT[i]);
	free(yx);
	
	return;
}

void EvaluateReservoirTrain(struct Reservoir *rc, char *filename) {
	int i, j, sum_node = rc->node_num * 2;
	double NMSE, sum_y, sqsum_y, variance_y, one_y, node_vc_difference;
	struct Histdata hist;
	FILE *fp_out;
	
	for(i = 0, NMSE = sum_y = sqsum_y = 0.0; i < rc->trainingNum; i++) {
		sum_y += rc->target_train[i];
		sqsum_y += rc->target_train[i] * rc->target_train[i];
		for(j = 0, one_y = 0.0; j < sum_node; j++) one_y += rc->weight[j] * rc->train_node[i][j];
		NMSE += (rc->target_train[i] - one_y) * (rc->target_train[i] - one_y);
		rc->output_train[i] = one_y;
	}
	
	variance_y = (sqsum_y - sum_y * sum_y / (double)rc->trainingNum) / (double)rc->trainingNum;
	NMSE /= variance_y * (double)rc->trainingNum;
	printf("Training result: Variance = %e, NMSE = %e\n", variance_y, NMSE);
	
	if(filename[0] != '0') {
		fp_out = fopen(filename, "w");
		fprintf(fp_out, "Index\tInput Santa Fe\tTarget\tTraining result\tNode vc difference\tVar of nodes at an each data\n");
		for(i = 0; i < rc->trainingNum; i++) {
			SetHistgramParameter(&hist, -1.0, 2.0, 100);
			for(j = 0, node_vc_difference = 0.0; j < rc->node_num; j++) {
				node_vc_difference +=  (rc->train_node[0][j] - rc->train_node[i][j]) * (rc->train_node[0][j] - rc->train_node[i][j]);
				node_vc_difference +=  (rc->train_node[0][rc->node_num + j] - rc->train_node[i][rc->node_num + j]) * (rc->train_node[0][rc->node_num + j] - rc->train_node[i][rc->node_num + j]);
				InputHistgram(&hist, rc->train_node[i][j]);
				InputHistgram(&hist, rc->train_node[i][rc->node_num + j]);
			}
			OutputHistgram(&hist, "0");
			fprintf(fp_out, "%d\t%e\t%e\t%e\t%e\n", i, rc->input_train[i], rc->target_train[i], rc->output_train[i], sqrt(node_vc_difference), hist.variance);
		}
		fclose(fp_out);
	}
	
	rc->NMSE_train = NMSE;
	
	return;
}

void EvaluateReservoirTest(struct Reservoir *rc, char *filename) {
	int i, j, sum_node = rc->node_num * 2;
	double NMSE, sum_y, sqsum_y, variance_y, one_y, node_vc_difference;
	struct Histdata hist;
	FILE *fp_out;
	
	for(i = 0, NMSE = sum_y = sqsum_y = 0.0; i < rc->testNum; i++) {
		sum_y += rc->target_test[i];
		sqsum_y += rc->target_test[i] * rc->target_test[i];
		for(j = 0, one_y = 0.0; j < sum_node; j++) one_y += rc->weight[j] * rc->test_node[i][j];
		NMSE += (rc->target_test[i] - one_y) * (rc->target_test[i] - one_y);
		rc->output_test[i] = one_y;
	}
	
	variance_y = (sqsum_y - sum_y * sum_y / (double)rc->testNum) / (double)rc->testNum;
	NMSE /= variance_y * (double)rc->testNum;
	printf("Test result: Variance = %e, NMSE = %e\n", variance_y, NMSE);
	
	if(filename[0] != '0') {
		fp_out = fopen(filename, "w");
		fprintf(fp_out, "Index\tInput Santa Fe\tTarget\tTest result\tNode vc difference\tVar of nodes at an each data\n");
		for(i = 0; i < rc->testNum; i++) {
			SetHistgramParameter(&hist, -1.0, 2.0, 100);
			for(j = 0, node_vc_difference = 0.0; j < rc->node_num; j++) {
				node_vc_difference +=  (rc->test_node[0][j] - rc->test_node[i][j]) * (rc->test_node[0][j] - rc->test_node[i][j]);
				node_vc_difference +=  (rc->test_node[0][rc->node_num + j] - rc->test_node[i][rc->node_num + j]) * (rc->test_node[0][rc->node_num + j] - rc->test_node[i][rc->node_num + j]);
				InputHistgram(&hist, rc->test_node[i][j]);
				InputHistgram(&hist, rc->test_node[i][rc->node_num + j]);
			}
			OutputHistgram(&hist, "0");
			fprintf(fp_out, "%d\t%e\t%e\t%e\t%e\t%e\n", i, rc->input_test[i], rc->target_test[i], rc->output_test[i], sqrt(node_vc_difference), hist.variance);
		}
		fclose(fp_out);
	}
	
	rc->NMSE_test = NMSE;
	
	return;
}

void Surrogate_Node(struct Reservoir *rc) {
	int i, j, rand_num[rc->node_num];
	double tmp;
	
	init_gen_rand(333);
	for(i = 0; i < rc->node_num; i++) rand_num[i] = (int)(genrand_open_open() * rc->node_num);
	
	for(i = 0; i < rc->trainingNum; i++) {
		for(j = 0; j < rc->node_num; j++) {
			tmp = rc->train_node[i][j];
			rc->train_node[i][j] = rc->train_node[i][rand_num[j]];
			rc->train_node[i][rand_num[j]] = tmp;
		}
	}
	
	return;
}

void CalcLaserEOfeedback(struct LaserEOfeedback *laser, double input_MZM) {
	int i, j, k;
	double a[2][2], ini[2][2], b[2][2][4], cosTheta, feedback[2], interval_L = h / laser[0].p.tau_L, interval_H = h / laser[0].p.tau_H;
	
	for(i = 0; i < 2; i++) {
		ini[i][0] = laser[i].x;
		ini[i][1] = laser[i].y;
		if(i == 0) cosTheta = cos(laser[0].p.kappa * laser[0].p.sigma * laser[0].xdelay[laser[0].delay_index[0]] + (1.0 - laser[0].p.kappa) * laser[1].p.sigma * laser[1].xdelay[laser[1].delay_index[1]] + laser[0].p.phi);
		else cosTheta = cos(laser[1].p.kappa * laser[1].p.sigma * laser[1].xdelay[laser[1].delay_index[0]] + (1.0 - laser[1].p.kappa) * laser[0].p.sigma * laser[0].xdelay[laser[0].delay_index[1]] + laser[1].p.phi);
		//feedback[i] = laser[i].p.beta * cosTheta * cosTheta;
		feedback[i] = laser[i].p.beta * 2.0 * input_MZM * cosTheta * cosTheta;// 追加分
	}
	
	for(i = 0; i < 4; i++) {
		if(i == 0) for(j = 0; j < 2; j++) for(k = 0; k < 2; k++) a[j][k] = ini[j][k];
		else if(i == 1) for(j = 0; j < 2; j++) for(k = 0; k < 2; k++) a[j][k] = ini[j][k] + b[j][k][0] / 2.0;
		else if(i == 2) for(j = 0; j < 2; j++) for(k = 0; k < 2; k++) a[j][k] = ini[j][k] + b[j][k][1] / 2.0;
		else if(i == 3) for(j = 0; j < 2; j++) for(k = 0; k < 2; k++) a[j][k] = ini[j][k] + b[j][k][2];
		
		for(j = 0; j < 2; j++) {
			b[j][0][i] = interval_L * ( - (1.0 + laser[j].p.epsilon) * a[j][0] - a[j][1] + feedback[j]);
			b[j][1][i] = interval_H * a[j][0];
		}
	}
	
	for(i = 0; i < 2; i++) for(j = 0; j < 2; j++) a[i][j] = ini[i][j] + (b[i][j][0] + 2.0 * b[i][j][1] + 2.0 * b[i][j][2] + b[i][j][3]) / 6.0;
	
	for(i = 0; i < 2; i++) {
		laser[i].xdelay[laser[i].delay_index[laser[i].long_delay]] = laser[i].x;
		for(j = 0; j < 2; j++) laser[i].delay_index[j] = (laser[i].delay_index[j] + 1) % laser[i].descretized_num[laser[i].long_delay];
		a[i][0] += sqrt( - 2.0 * WGSD * h * log(genrand_open_open()) ) * cos(2.0 * M_PI * genrand_open_open());
		laser[i].x = a[i][0];
		laser[i].y = a[i][1];
	}
	
	return;
}

void SetHistgramParameter(struct Histdata *histgram, double min, double max, int bin) {
	int i;
	
	histgram->data = (int*)malloc(sizeof(int) * bin);
	if(histgram->data == NULL) fprintf(stderr, "error : histgram\n");
	histgram->min = min;
	histgram->max = max;
	histgram->bin = bin;
	histgram->div = (max - min) / bin;
	for(i = 0; i < bin; i++) histgram->data[i] = 0;
	histgram->minData = histgram->max;
	histgram->maxData = histgram->min;
	histgram->average = 0.0;
	histgram->variance = 0.0;
	histgram->n = 0;
	
	return;
}

void InputHistgram(struct Histdata *histgram, double inputData) {
	
	if(inputData >= histgram->min && inputData < histgram->max) histgram->data[(int)((inputData - histgram->min) / histgram->div + EPS)]++;
	if(inputData < histgram->minData) histgram->minData = inputData;
	if(inputData > histgram->maxData) histgram->maxData = inputData;
	histgram->average += inputData;
	histgram->variance += inputData * inputData;
	histgram->n++;
	
	return;
}

void OutputHistgram(struct Histdata *histgram, char *filename) {
	int i, maxDatanum;
	FILE *outputfp;
	
	histgram->variance = (histgram->variance - histgram->average * histgram->average / histgram->n) / histgram->n;
	histgram->average /= histgram->n;
	for(i = 0, maxDatanum = 0, histgram->entropy = 0.0; i < histgram->bin; i++) {
		if(histgram->data[i] > maxDatanum) {
			maxDatanum = histgram->data[i];
			histgram->median = (i + 0.5) * histgram->div + histgram->min;
		}
		if(histgram->data[i] != 0.0) histgram->entropy -= (double)histgram->data[i] / (double)(histgram->n) * log2((double)histgram->data[i] / (double)(histgram->n));
	}
	histgram->entropy /= -log2(1.0 / (double)histgram->bin);
	
	if(filename[0] != '0') {
		outputfp = fopen(filename, "w");
		fprintf(outputfp, "Maximum data : %e\nMinimum data : %e\nMedian : %e\nAverage : %e\nVariance : %e\nNormalized entropy : %e\n", histgram->maxData, histgram->minData, histgram->median, histgram->average, histgram->variance, histgram->entropy);
		fprintf(outputfp, "Data\tProbability\tThe number of points\n");
		for(i = 0; i < histgram->bin; i++) fprintf(outputfp, "%e\t%e\t%d\n", (i + 0.5) * histgram->div + histgram->min, (double)histgram->data[i] / (double)(histgram->n), histgram->data[i]);
		fclose(outputfp);
	}
	
	free(histgram->data);
	
	return;
}

void Set_CC_fixed_shift(struct cc_data *cc) {
	
	cc->sum[0] = cc->sum[1] = 0.0;
	cc->sq_sum[0] = cc->sq_sum[1] = 0.0;
	cc->product_sum = 0.0;
	cc->difference = 0.0;
	cc->num = 0;
	
	return;
}

void Input_CC_fixed_shift(struct cc_data *cc, double input1, double input2) {
	
	cc->sum[0] += input1;
	cc->sum[1] += input2;
	cc->sq_sum[0] += input1 * input1;
	cc->sq_sum[1] += input2 * input2;
	cc->product_sum += input1 * input2;
	cc->difference += fabs(input1 - input2);
	cc->num++;
	
	return;
}

void Output_CC_fixed_shift(struct cc_data *cc) {
	double var[2];
	int i;
	
	for(i = 0; i < 2; i++) var[i] = (cc->sq_sum[i] - cc->sum[i] * cc->sum[i] / (double)cc->num) / (double)cc->num;
	cc->correlation = (cc->product_sum - cc->sum[1] / (double)cc->num * cc->sum[0]) / (sqrt(var[0] * var[1]) * (double)cc->num);
	cc->difference /= (double)cc->num;
	
	return;
}
