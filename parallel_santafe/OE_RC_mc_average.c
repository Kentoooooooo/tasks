#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <windows.h>
#include "dSFMT.h"

#define h 1e-10
#define dh 1e-4
#define EPS 1e-8
#define WGSD 1.0e4

struct parameter {
	double tau_L, tau_H, delay[2], beta, kappa, sigma, phi, epsilon;
};

struct LaserEOfeedback {
	struct parameter p;
	double x, y, *xdelay;
	int descretized_num[2], delay_index[2], long_delay;
};

struct Reservoir {
	int trainingNum, testNum, mask_num, node_num, allNum, mask_interval, node_interval, input_interval, maskIndex, nodeIndex, inputIndex, nodeIndex_l2, inputIndex_l2;
	double gamma, *mask, **train_node, **test_node, *target_train, *target_test, *weight, *input_train, *input_test, *output_train, *output_test, input_signal, input_MZM, input_MZM_phi;
	double NMSE_train, NMSE_test;
};

struct Histdata {
	double min, max, div, minData, maxData, median, average, variance, entropy;
	int bin, *data, n;
};

void CalcLaserEOfeedback(struct LaserEOfeedback *laser, double input_signal);
void SetRC_Param(struct Reservoir *rc, double mask_interval, int mask_num, int mask_node_ratio);
double *MakeMask(int mask_num, int seed_mask);
void ReadSantaFe(struct Reservoir *rc);
void QR(struct Reservoir *rc);
void EvaluateReservoirTrain(struct Reservoir *rc);
void EvaluateReservoirTest(struct Reservoir *rc);
void Surrogate_Node(struct Reservoir *rc);
void free_rc_memory(struct Reservoir *rc);
void SetHistgramParam(struct Histdata *histgram, double min, double max, int bin);
void InputHistgram(struct Histdata *histgram, double inputData);
void OutputHistgram(struct Histdata *histgram, char *filename);

int main(int argc, char *argv[]) {
	struct LaserEOfeedback laser[2];
	struct Reservoir rc;
	struct Histdata hist[6], hist_NMSE_train, hist_NMSE_test, hist_SNR_laser1, hist_SNR_laser2, hist_node_state_sd_laser1, hist_node_state_sd_laser2, hist_node_state_sd_all, hist_node_state_entropy[3];
	double param, param_min, param_max, param_interval, t_transient, t_plot;
	double buff[8], *wfmData[8];
	int i, j, k, n, p, s, t, transient, calc_num, divi, out, para_flag, iteration, plotnum, feedback_memory, index, ave_idx, ave_num;
	char buffer_ch[128], filename[128];
	clock_t start, end;
	FILE *inputfp, *outputfp[8];
	
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
		fscanf(inputfp, "%d", &ave_num);
		
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
		if(para_flag == 11) sprintf(buffer_ch, "Delay time (MZM2) and propagation time (MZM1) \\ft\\n\\d2\\n [ms]");
		if(para_flag == 12) sprintf(buffer_ch, "Mask pattern");
		outputfp[0] = fopen("param_change.txt", "w");
		fprintf(outputfp[0], "%s\tNMSE for training\tVariance of NMSE for training\tNMSE for test\tVariance of NMSE for test\tSNR for laser1\tVariance of SNR for laser1\tSNR for laser2\tVariance of SNR for laser2\tAve. of var. of node states\tAve. of node state entropy\n", buffer_ch);
		outputfp[1] = fopen("param_change_laser1_node_state_sd.txt", "w");
		outputfp[2] = fopen("param_change_laser2_node_state_sd.txt", "w");
		outputfp[3] = fopen("param_change_all_node_state_sd.txt", "w");
		outputfp[4] = fopen("param_change_laser1_node_state_entropy.txt", "w");
		outputfp[5] = fopen("param_change_laser2_node_state_entropy.txt", "w");
		outputfp[6] = fopen("param_change_all_node_state_entropy.txt", "w");
		for(i = 1; i < 7; i++) {
			fprintf(outputfp[i], "%s", buffer_ch);
			for(j = 0; j < ave_num; j++) fprintf(outputfp[i], "\t%d", j);
		}
		for(i = 1; i < 4; i++) fprintf(outputfp[i], "\tAve. of node states sd\tVar. of sd\n");
		for(i = 4; i < 7; i++) fprintf(outputfp[i], "\tAve. of node states entropy\tVar. of entropy\n");
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
	
	SetRC_Param(&rc, 0.2e-6, 50, 1);
	
	calc_num = rc.allNum * rc.input_interval;
	if(laser[0].p.delay[0] > laser[1].p.delay[0]) calc_num++;
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
				calc_num = rc.allNum * rc.input_interval;
				if(laser[0].p.delay[0] > laser[1].p.delay[0]) calc_num++;
				printf("\nNum of nodes = %d (node interval = %e)\n", rc.node_num, rc.node_interval * h);
			}
			if(para_flag == 8) {
				free_rc_memory(&rc);
				SetRC_Param(&rc, 0.4e-6, 25, 2);
				calc_num = rc.allNum * rc.input_interval;
				if(laser[0].p.delay[0] > laser[1].p.delay[0]) calc_num++;
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
		
		SetHistgramParam(&hist_NMSE_train, 0.0, 1.0, 100);
		SetHistgramParam(&hist_NMSE_test, 0.0, 1.0, 100);
		SetHistgramParam(&hist_SNR_laser1, 0.0, 40.0, 100);
		SetHistgramParam(&hist_SNR_laser2, 0.0, 40.0, 100);
		SetHistgramParam(&hist_node_state_sd_laser1, 0.0, 1.0, 100);
		SetHistgramParam(&hist_node_state_sd_laser2, 0.0, 1.0, 100);
		SetHistgramParam(&hist_node_state_sd_all, 0.0, 1.0, 100);
		for(i = 0; i < 3; i++) SetHistgramParam(&hist_node_state_entropy[i], 0.0, 1.0, 100);
		
		if(out == 2) for(i = 1; i < 7; i++) fprintf(outputfp[i], "%e", param);
		
		for(ave_idx = 0; ave_idx < ave_num; ave_idx++) {
			
			if(rc.mask != NULL) free(rc.mask);
			rc.mask = MakeMask(rc.mask_num, 1234 + ave_idx);
			
			for(i = 0; i < 2; i++) {
				laser[i].x = 0.05 + (double)i * 0.05;
				laser[i].y = 0.0;
				for(j = 0; j < laser[i].descretized_num[laser[i].long_delay]; j++) laser[i].xdelay[j] = 0.05 + (double)i * 0.05;
				for(j = 0; j < 2; j++) laser[i].delay_index[j] = laser[i].descretized_num[laser[i].long_delay] - laser[i].descretized_num[j];
			}
			init_gen_rand(222 + ave_idx);
			
			SetHistgramParam(&hist[0], -0.7, 1.5, 220);
			SetHistgramParam(&hist[1], -0.7, 1.5, 220);
			SetHistgramParam(&hist[2], -0.7, 1.5, 220);
			SetHistgramParam(&hist[3], -0.7, 1.5, 220);
			SetHistgramParam(&hist[6], -0.7, 1.5, 220);
			
			for(i = 0; i < transient; i++) CalcLaserEOfeedback(laser, 0.5);
			
			calc_num = rc.trainingNum * rc.input_interval;
			rc.maskIndex = rc.nodeIndex = rc.inputIndex = 0;
			
			for(i = 1; i <= calc_num; i++) {
				
				rc.input_signal = rc.gamma * rc.mask[rc.maskIndex] * rc.input_train[rc.inputIndex];
				rc.input_MZM = cos(rc.input_signal * M_PI / 4.0 + rc.input_MZM_phi);
				rc.input_MZM *= rc.input_MZM;
				
				// CalcLaserEOfeedback(laser, 0.5);
				CalcLaserEOfeedback(laser, rc.input_MZM);
				InputHistgram(&hist[0], laser[0].x);
				InputHistgram(&hist[1], laser[1].x);
				
				if(i % rc.node_interval == rc.node_interval / 2) {
					rc.train_node[rc.inputIndex][rc.nodeIndex] = laser[0].x;
					rc.train_node[rc.inputIndex][rc.node_num + rc.nodeIndex] = laser[1].x;
					rc.nodeIndex = (rc.nodeIndex + 1) % rc.node_num;
					InputHistgram(&hist[2], laser[0].x);
					InputHistgram(&hist[3], laser[1].x);
					InputHistgram(&hist[6], laser[0].x);
					InputHistgram(&hist[6], laser[1].x);
				}
				
				if(i % rc.mask_interval == 0) rc.maskIndex = (rc.maskIndex + 1) % rc.mask_num;
				if(i % rc.input_interval == 0) rc.inputIndex++;
			}
			
			calc_num = rc.testNum * rc.input_interval;
			rc.maskIndex = rc.nodeIndex = rc.inputIndex = 0;
			
			for(i = 1; i <= calc_num; i++) {
				
				rc.input_signal = rc.gamma * rc.mask[rc.maskIndex] * rc.input_test[rc.inputIndex];
				rc.input_MZM = cos(rc.input_signal * M_PI / 4.0 + rc.input_MZM_phi);
				rc.input_MZM *= rc.input_MZM;
				
				// CalcLaserEOfeedback(laser, 0.5);
				CalcLaserEOfeedback(laser, rc.input_MZM);
				InputHistgram(&hist[0], laser[0].x);
				InputHistgram(&hist[1], laser[1].x);
				
				if(i % rc.node_interval == rc.node_interval / 2) {
					rc.test_node[rc.inputIndex][rc.nodeIndex] = laser[0].x;
					rc.test_node[rc.inputIndex][rc.node_num + rc.nodeIndex] = laser[1].x;
					rc.nodeIndex = (rc.nodeIndex + 1) % rc.node_num;
					InputHistgram(&hist[2], laser[0].x);
					InputHistgram(&hist[3], laser[1].x);
					InputHistgram(&hist[6], laser[0].x);
					InputHistgram(&hist[6], laser[1].x);
				}
				
				if(i % rc.mask_interval == 0) rc.maskIndex = (rc.maskIndex + 1) % rc.mask_num;
				if(i % rc.input_interval == 0) rc.inputIndex++;
			}
			
			SetHistgramParam(&hist[4], -0.7, 1.5, 220);
			SetHistgramParam(&hist[5], -0.7, 1.5, 220);
			for(i = 0; i < transient; i++) CalcLaserEOfeedback(laser, 0.5);
			for(i = 0; i < transient; i++) {
				CalcLaserEOfeedback(laser, 0.5);
				InputHistgram(&hist[4], laser[0].x);
				InputHistgram(&hist[5], laser[1].x);
			}
			
			OutputHistgram(&hist[0], "0");
			OutputHistgram(&hist[1], "0");
			OutputHistgram(&hist[2], "0");
			OutputHistgram(&hist[3], "0");
			OutputHistgram(&hist[4], "0");
			OutputHistgram(&hist[5], "0");
			OutputHistgram(&hist[6], "0");
			printf("SNR for laser1 = %2.3f\n", 10.0 * log10(hist[0].variance / hist[4].variance));
			printf("SNR for laser2 = %2.3f\n", 10.0 * log10(hist[1].variance / hist[5].variance));
			
			// Surrogate_Node(&rc);
			QR(&rc);
			EvaluateReservoirTrain(&rc);
			EvaluateReservoirTest(&rc);
			printf("NMSE for training = %1.4f\tNMSE for test = %1.4f\n", rc.NMSE_train, rc.NMSE_test);
			
			InputHistgram(&hist_NMSE_train, rc.NMSE_train);
			InputHistgram(&hist_NMSE_test, rc.NMSE_test);
			InputHistgram(&hist_SNR_laser1, 10.0 * log10(hist[0].variance / hist[4].variance));
			InputHistgram(&hist_SNR_laser2, 10.0 * log10(hist[1].variance / hist[5].variance));
			InputHistgram(&hist_node_state_sd_laser1, sqrt(hist[2].variance));
			InputHistgram(&hist_node_state_sd_laser2, sqrt(hist[3].variance));
			InputHistgram(&hist_node_state_sd_all, sqrt(hist[6].variance));
			InputHistgram(&hist_node_state_entropy[0], hist[2].entropy);
			InputHistgram(&hist_node_state_entropy[1], hist[3].entropy);
			InputHistgram(&hist_node_state_entropy[2], hist[6].entropy);
			
			if(out == 2) {
				fprintf(outputfp[1], "\t%e", sqrt(hist[2].variance));
				fprintf(outputfp[2], "\t%e", sqrt(hist[3].variance));
				fprintf(outputfp[3], "\t%e", sqrt(hist[6].variance));
				fprintf(outputfp[4], "\t%e", hist[2].entropy);
				fprintf(outputfp[5], "\t%e", hist[3].entropy);
				fprintf(outputfp[6], "\t%e", hist[6].entropy);
			}
			
			end = clock();
			printf("%f\n", (double)(end - start) / CLOCKS_PER_SEC);
		}
		
		OutputHistgram(&hist_NMSE_train, "0");
		OutputHistgram(&hist_NMSE_test, "0");
		OutputHistgram(&hist_SNR_laser1, "0");
		OutputHistgram(&hist_SNR_laser2, "0");
		OutputHistgram(&hist_node_state_sd_laser1, "0");
		OutputHistgram(&hist_node_state_sd_laser2, "0");
		OutputHistgram(&hist_node_state_sd_all, "0");
		for(i = 0; i < 3; i++) OutputHistgram(&hist_node_state_entropy[i], "0");
		
		if(out == 1) {
			if(para_flag == 2) {
				
				outputfp[0] = fopen("rc_training_result.txt", "w");
				fprintf(outputfp[0], "Index\tInput Santa Fe\tTarget\tTraining result\n");
				for(i = 0; i < rc.trainingNum; i++) fprintf(outputfp[0], "%d\t%e\t%e\t%e\n", i, rc.input_train[i], rc.target_train[i], rc.output_train[i]);
				fclose(outputfp[0]);
				
				outputfp[0] = fopen("rc_test_result.txt", "w");
				fprintf(outputfp[0], "Index\tInput Santa Fe\tTarget\tTest result\n");
				for(i = 0; i < rc.testNum; i++) fprintf(outputfp[0], "%d\t%e\t%e\t%e\n", i, rc.input_test[i], rc.target_test[i], rc.output_test[i]);
				fclose(outputfp[0]);
			}
			break;
		}
		if(out == 2) {
			fprintf(outputfp[0], "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", param, hist_NMSE_train.average, hist_NMSE_train.variance, hist_NMSE_test.average, hist_NMSE_test.variance, hist_SNR_laser1.average, hist_SNR_laser1.variance, hist_SNR_laser2.average, hist_SNR_laser2.variance, hist_node_state_sd_all.average, hist_node_state_entropy[2].average);
			fprintf(outputfp[1], "\t%e\t%e\n", hist_node_state_sd_laser1.average, hist_node_state_sd_laser1.variance);
			fprintf(outputfp[2], "\t%e\t%e\n", hist_node_state_sd_laser2.average, hist_node_state_sd_laser2.variance);
			fprintf(outputfp[3], "\t%e\t%e\n", hist_node_state_sd_all.average, hist_node_state_sd_all.variance);
			fprintf(outputfp[4], "\t%e\t%e\n", hist_node_state_entropy[0].average, hist_node_state_entropy[0].variance);
			fprintf(outputfp[5], "\t%e\t%e\n", hist_node_state_entropy[1].average, hist_node_state_entropy[1].variance);
			fprintf(outputfp[6], "\t%e\t%e\n", hist_node_state_entropy[2].average, hist_node_state_entropy[2].variance);
		}
	}
	
	free_rc_memory(&rc);
	
	for(i = 0; i < 2; i++) free(laser[i].xdelay);
	if(out == 2) for(i = 0; i < 7; i++) fclose(outputfp[i]);
	
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
	rc->gamma = 1.0 / 255.0;
	rc->input_MZM_phi = - M_PI / 4.0;
	
	// SantaFe
	rc->target_train = (double*)malloc(sizeof(double) * rc->trainingNum);
	rc->target_test = (double*)malloc(sizeof(double) * rc->testNum);
	rc->input_train = (double*)malloc(sizeof(double) * rc->trainingNum);
	rc->input_test = (double*)malloc(sizeof(double) * rc->testNum);
	rc->output_train = (double*)malloc(sizeof(double) * rc->trainingNum);
	rc->output_test = (double*)malloc(sizeof(double) * rc->testNum);
	ReadSantaFe(rc);
	
	rc->mask_num = mask_num;
	rc->mask_interval = (int)(mask_interval / h + EPS);
	rc->input_interval = rc->mask_interval * rc->mask_num;
	
	rc->node_num = rc->mask_num * mask_node_ratio;
	rc->node_interval = rc->mask_interval / mask_node_ratio;
	
	// Make a mask for reservoir computing
	rc->mask = MakeMask(rc->mask_num, 1234);
	
	rc->train_node = (double**)malloc(sizeof(double*) * rc->trainingNum);
	for(i = 0; i <  rc->trainingNum; i++) rc->train_node[i] = (double*)malloc(sizeof(double) * rc->node_num * 2);
	rc->test_node = (double**)malloc(sizeof(double*) * rc->testNum);
	for(i = 0; i <  rc->testNum; i++) rc->test_node[i] = (double*)malloc(sizeof(double) * rc->node_num * 2);
	rc->weight = (double*)malloc(sizeof(double) * rc->node_num * 2);
	
	return;
}

double *MakeMask(int mask_num, int seed_mask) {
	double *mask;
	int i, rand_num, old_rand_num, mask_pattern_value;
	
	mask_pattern_value = 4.0;
	old_rand_num = -1;
	
	mask = (double*)malloc(sizeof(double) * mask_num);
	
	init_gen_rand(seed_mask);
	for(i = 0; i < mask_num; i++) {
		do {
			rand_num = (int)(mask_pattern_value * genrand_open_open());
		} while(rand_num == old_rand_num);
		old_rand_num = rand_num;
		
		if(rand_num == 0) mask[i] = -1.0;
		if(rand_num == 1) mask[i] = -0.3;
		if(rand_num == 2) mask[i] = 0.3;
		if(rand_num == 3) mask[i] = 1.0;
	}
	
	return mask;
}

void ReadSantaFe(struct Reservoir *rc) {
	int i, count, ret, data, tmp[5000];
	FILE *inputfp;
	
	inputfp = fopen("santafe.txt", "r");
	if(inputfp == NULL) {
		fprintf(stderr, "Can not open the file of santafe.txt.");
		exit(2);
	}
	
	count = 0;
	while((ret = fscanf(inputfp, "%d", &data)) != EOF) {
		tmp[count] = (double)data;
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

void EvaluateReservoirTrain(struct Reservoir *rc) {
	int i, j, sum_node = rc->node_num * 2;
	double NMSE, sum_y, sqsum_y, variance_y, one_y;
	
	for(i = 0, NMSE = sum_y = sqsum_y = 0.0; i < rc->trainingNum; i++) {
		sum_y += rc->target_train[i];
		sqsum_y += rc->target_train[i] * rc->target_train[i];
		for(j = 0, one_y = 0.0; j < sum_node; j++) one_y += rc->weight[j] * rc->train_node[i][j];
		NMSE += (rc->target_train[i] - one_y) * (rc->target_train[i] - one_y);
		rc->output_train[i] = one_y;
	}
	
	variance_y = (sqsum_y - sum_y * sum_y / (double)rc->trainingNum) / (double)rc->trainingNum;
	NMSE /= variance_y * (double)rc->trainingNum;
	// printf("Training result: Variance = %e, NMSE = %e\n", variance_y, NMSE);
	rc->NMSE_train = NMSE;
	
	return;
}

void EvaluateReservoirTest(struct Reservoir *rc) {
	int i, j, sum_node = rc->node_num * 2;
	double NMSE, sum_y, sqsum_y, variance_y, one_y;
	
	for(i = 0, NMSE = sum_y = sqsum_y = 0.0; i < rc->testNum; i++) {
		sum_y += rc->target_test[i];
		sqsum_y += rc->target_test[i] * rc->target_test[i];
		for(j = 0, one_y = 0.0; j < sum_node; j++) one_y += rc->weight[j] * rc->test_node[i][j];
		NMSE += (rc->target_test[i] - one_y) * (rc->target_test[i] - one_y);
		rc->output_test[i] = one_y;
	}
	
	variance_y = (sqsum_y - sum_y * sum_y / (double)rc->testNum) / (double)rc->testNum;
	NMSE /= variance_y * (double)rc->testNum;
	// printf("Test result: Variance = %e, NMSE = %e\n", variance_y, NMSE);
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

void SetHistgramParam(struct Histdata *histgram, double min, double max, int bin) {
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
