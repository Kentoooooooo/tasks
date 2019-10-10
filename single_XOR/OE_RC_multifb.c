#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <windows.h>
#include "dSFMT.h"

#define h 1e-10 //stride
#define dh 1e-4
#define EPS 1e-8
#define WGSD 1.0e4
#define OUT_NODE_STATE 0

struct parameter {
	double tau_L, tau_H, delay[10], beta, phi, epsilon, weight_fb[10], fb_str;
};

struct LaserEOfeedback {
	struct parameter p;
	double x, y, *xdelay;
	int descretized_num[10], fbNum, fbIndex[10], long_fb;
};

struct Reservoir {
	int trainingNum, testNum, mask_num, node_num, allNum, mask_interval, node_interval, input_interval, maskIndex, nodeIndex, inputIndex, delay_step;
	double gamma, *mask, **nodeState, *target, *weight, *inputNARMA10, *inputdata, *output, input_signal;
};

struct Histdata {
	double min, max, div, minData, maxData, median, average, var, entropy;
	int bin, *data, n;
};

void CalcLaserEOfeedback(struct LaserEOfeedback *laser, double input_signal);
void SetRC_Param(struct Reservoir *rc, double mask_interval, int mask_num, int mask_node_ratio);
void MakeMask(struct Reservoir *rc);
double calc_ser(double *target, double *signal, int num);
void ReadSantaFe(struct Reservoir *rc);
void QR(struct Reservoir *rc);
double EvaluateReservoirTrain(struct Reservoir *rc);
double EvaluateReservoirTest(struct Reservoir *rc);
void free_rc_memory(struct Reservoir *rc);
void SetHistgramParameter(struct Histdata *histgram, double min, double max, int bin);
void InputHistgram(struct Histdata *histgram, double inputData);
void OutputHistgram(struct Histdata *histgram, char *filename);

int main(int argc, char *argv[]) {
	struct LaserEOfeedback laser;
	struct Reservoir rc;
	struct Histdata hist;
	double param, param_min, param_max, param_interval, max_delay, t_transient, t_plot;
	double buff[8], *wfmData[8];
	int i, j, k, n, p, s, t, transient, calc_num, divi, out, para_flag, iteration, plotnum, feedback_memory;
	char buffer_ch[128], filename[128];
	clock_t start, end;
	FILE *inputfp, *outputfp[4];

	//	input data from file(argv[1])
	if (argc == 2) {
		inputfp = fopen(argv[1], "r");
		if(inputfp == NULL) {
			fprintf(stderr, "%sを開くことができません", argv[1]);
			exit(2);
		}
		fscanf(inputfp, "%lf", &(laser.p.tau_L));
		fscanf(inputfp, "%lf", &(laser.p.tau_H));
		fscanf(inputfp, "%d", &(laser.fbNum));
		for(i = 0; i < laser.fbNum; i++) fscanf(inputfp, "%lf", &(laser.p.delay[i]));
		fscanf(inputfp, "%lf", &(laser.p.beta));
		fscanf(inputfp, "%lf", &(laser.p.fb_str));
		fscanf(inputfp, "%lf", &(laser.p.phi));
		fscanf(inputfp, "%lf %lf %lf", &param_min, &param_max, &param_interval);
		fscanf(inputfp, "%lf %lf", &t_transient, &t_plot);
		fscanf(inputfp, "%d", &divi);
		fscanf(inputfp, "%d", &(rc.delay_step));
		fscanf(inputfp, "%lf", &(rc.gamma));

		fclose(inputfp);
	} else {
		fprintf(stderr, "引数のエラーです\n");
		exit(1);
	}

	const int amp_max_int = (param_max * param_interval >= 0.0) ? (int)(param_max / param_interval + EPS) : (int)(param_max / param_interval - EPS);
	const int amp_min_int = (param_min * param_interval >= 0.0) ? (int)(param_min / param_interval + EPS) : (int)(param_min / param_interval - EPS);

	//	selection of calculation mode
	printf("Select a calculation-mode (1 or 2)\n");
	printf("\t 1 --> Temporal waveforms\n");
	printf("\t 2 --> Parameter change\n");
	printf("inputfp --> ");
	scanf("%d", &out);
	if(out == 1) {
		printf("Select a output-mode (1 or 2 or 3)\n");
		printf("\t 1 --> Temporal waveforms of x(t) and y(t)\n");
		printf("\t 2 --> Reservoir Computing\n");
		printf("inputfp --> ");
		scanf("%d", &para_flag);
		if(para_flag == 1) {
			outputfp[0] = fopen("twfm.txt", "w");
			fprintf(outputfp[0], "Time [\\fm\\ns]\tx(t)\tMask\tIdeal\tInput\tNode\n");
		} else if(para_flag == 2) {
		} else {
			fprintf(stderr, "error : para_flag = %d\n", para_flag);
			exit(1);
		}
	//selection of chg_parameter
	} else if(out == 2) {
		printf("Select a parameter which you change (1 or 2 or 3)\n");
		printf("\t 1 --> Beta\n");
		printf("\t 2 --> Feedback delay time 1\n");
		printf("\t 3 --> Bias point for MZM\n");
		printf("\t 4 --> Gamma\n");
		printf("\t 5 --> Num of nodes\n");
		printf("\t 6 --> Feedback strength\n");
		printf("\t 7 --> Mask pattern\n");
		printf("\t 8 --> Nonlinearity Nu\n");
		printf("\t 9 --> Memory Tau_d\n");
		printf("inputfp --> ");
		scanf("%d", &para_flag);
		if(para_flag == 1) sprintf(buffer_ch, "Feedback strength \\fb\\n");
		if(para_flag == 2) sprintf(buffer_ch, "Feedback delay time 1 \\ft\\n\\d1\\n [\\fm\\ns]");
		if(para_flag == 3) sprintf(buffer_ch, "Bias point for MZM \\ff\\n");
		if(para_flag == 4) sprintf(buffer_ch, "Gamma \\fg\\n");
		if(para_flag == 5) sprintf(buffer_ch, "Number of nodes");
		if(para_flag == 6) sprintf(buffer_ch, "Feedback strength");
		if(para_flag == 7) sprintf(buffer_ch, "Mask pattern");
		if(para_flag == 8) sprintf(buffer_ch, "Delay step \\ft\\n\\dD\\n");

		outputfp[0] = fopen("parameter_change.txt", "w");
		fprintf(outputfp[0], "%s\tNMSE for training\tNMSE for test\tSER\tNormalized entropy\n", buffer_ch);

		if(OUT_NODE_STATE == 1) CreateDirectory(".\\node_state", NULL);
	} else {
		fprintf(stderr, "error : out = %d\n", out);
		exit(1);
	}

	/* Initialization... START */
	for(i = 0; i < laser.fbNum; i++) laser.p.weight_fb[i] = 1.0 / (double)(laser.fbNum);
	laser.p.phi *= M_PI;
	laser.p.tau_L = 1.0 / (2.0 * M_PI * laser.p.tau_L);
	laser.p.tau_H = 1.0 / (2.0 * M_PI * laser.p.tau_H);
	laser.p.epsilon = laser.p.tau_L / laser.p.tau_H;
	for(i = 0, laser.long_fb = 0; i < laser.fbNum; i++) {
		laser.descretized_num[i] = (int)(laser.p.delay[i] / h + EPS);
		if(laser.p.delay[i] > laser.p.delay[laser.long_fb]) laser.long_fb = i;
		printf("Delay time = %e (Descretized num = %d)\n", laser.p.delay[i], laser.descretized_num[i]);
	}
	laser.xdelay = (double*)malloc(sizeof(double) * laser.descretized_num[laser.long_fb]);
	printf("epsilon = %e, beta = %e\n", laser.p.epsilon, laser.p.beta);

	SetRC_Param(&rc, 0.2e-6, 50, 1);//mask_interval, mask_num, mask_node_ratio

	// SetRC_Param(&rc, 1.0e-6, 10, 10);
	// rc.mask[0] = -1.0;
	// rc.mask[1] = 1.0;
	// rc.mask[2] = -1.0;
	// rc.mask[3] = -1.0;
	// rc.mask[4] = 1.0;

	calc_num = rc.allNum * rc.input_interval;
	transient = (int)(t_transient / dh + EPS);
	plotnum = (int)(t_plot / dh + EPS);
	printf("Number of node          : %d\n", rc.node_num);
	printf("Number of one mask point: %d\n", rc.node_interval);
	printf("Bias point              : %e\n", laser.p.phi);
	printf("kappa                   : %e\n", laser.p.fb_str);
	printf("gamma                   : %e\n", rc.gamma);
	printf("delay_step              : %d\n", rc.delay_step);
	/* initialization... END */

	/* initialization for laser data */
	laser.x = 0.05;
	laser.y = 0.0;
	for(i = 0; i < laser.descretized_num[laser.long_fb]; i++) laser.xdelay[i] = 0.05;
	for(i = 0; i < laser.fbNum; i++) laser.fbIndex[i] = laser.descretized_num[laser.long_fb] - laser.descretized_num[i];
	/* initializatoin for laser data */

	//	MAIN calculation ----------------------
	for(s = amp_min_int; s <= amp_max_int; s++) {
		param = param_interval * (double)s;

		start = clock();

		if(out == 2) {
			if(para_flag == 1) {  //chg beta
				laser.p.beta = param;
				printf("\nbeta = %e\n", laser.p.beta);
			}
			if(para_flag == 2) {  //chg delay time
				laser.p.delay[0] = param;
				laser.descretized_num[0] = (int)(laser.p.delay[0] / h + EPS);
				for(i = 0, laser.long_fb = 0; i < laser.fbNum; i++) if(laser.p.delay[i] > laser.p.delay[laser.long_fb]) laser.long_fb = i;
				free(laser.xdelay);
				laser.xdelay = (double*)malloc(sizeof(double) * laser.descretized_num[laser.long_fb]);
				printf("\nDelay time 1 = %e (Descretized num = %d)\n", laser.p.delay[0], laser.descretized_num[0]);
				printf("Number of node: %d\n", rc.node_num);
				param *= 1.0e6;
			}
			if(para_flag == 3) {  //chg bias point
				laser.p.phi = M_PI * param;
				printf("\nBias point for MZM = %e * PI\n", param);
			}
			if(para_flag == 4) {  //chg gamma
				rc.gamma = param;
				printf("\ngamma = %e\n", rc.gamma);
			}
			if(para_flag == 5) { //chg num of nodes
				free_rc_memory(&rc);
				// SetRC_Param(&rc, 0.4e-6, 25, 2);
				SetRC_Param(&rc, 2.0e-6, 5, param / 5);
				rc.mask[0] = -0.1;
				rc.mask[1] = 0.1;
				rc.mask[2] = 0.1;
				rc.mask[3] = -0.1;
				rc.mask[4] = 0.1;
				calc_num = rc.allNum * rc.input_interval;
				printf("\nNum of nodes = %d (node interval = %e)\n", (int)param, rc.node_interval * h);
			}
			if(para_flag == 6) {  // chg fb_str
				laser.p.fb_str = param;
				printf("\nFeedback strength = %e\n", param);
			}
			if(para_flag == 7) {  //chg mask pattern
				printf("\nMask pattern = %e\n", param);
				for(i = 0, buff[0] = param; i < rc.mask_num; i++) {
					if((int)buff[0] % 2 == 1) rc.mask[i] = 0.1;
					else rc.mask[i] = -0.1;
					buff[0] = (int)(buff[0] / 2);
					printf("%+1.1f, ", rc.mask[i]);
				}
				printf("\n");
			}
			if(para_flag == 8) {  //chg delay_step
				rc.delay_step = param;

				SetRC_Param(&rc, 0.2e-6, 50, 1);

				printf("\ndelay_step = %d\n", (int)rc.delay_step);
			}
		}

		laser.x = 0.05;
		laser.y = 0.0;
		for(i = 0; i < laser.descretized_num[laser.long_fb]; i++) laser.xdelay[i] = 0.05;
		for(i = 0; i < laser.fbNum; i++) laser.fbIndex[i] = laser.descretized_num[laser.long_fb] - laser.descretized_num[i];
		init_gen_rand(222);

		SetHistgramParameter(&hist, -0.8, 0.8, 100);

		rc.maskIndex = rc.nodeIndex = rc.inputIndex = 0;

		iteration = -10;
		for(i = 0; i < 8; i++) buff[i] = 0.0;

		for(i = 0; i < transient; i++) CalcLaserEOfeedback(&laser, 0.0);

		for(i = 1; i <= 100 * rc.input_interval; i++) {

			rc.input_signal = rc.gamma * rc.mask[rc.maskIndex] * rc.inputdata[rc.inputIndex];

			CalcLaserEOfeedback(&laser, rc.input_signal);

			if(i % rc.mask_interval == 0) rc.maskIndex = (rc.maskIndex + 1) % rc.mask_num;
			if(i % rc.input_interval == 0) rc.inputIndex++;
		}

		for(i = 1; i <= calc_num; i++) {

			rc.input_signal = rc.gamma * rc.mask[rc.maskIndex] * rc.inputdata[rc.inputIndex];

			CalcLaserEOfeedback(&laser, rc.input_signal);

			if(out == 1) {
				if(para_flag == 1) {
					if(i % divi == 0) {
						if(i % rc.node_interval == rc.node_interval / 2)
							fprintf(outputfp[0], "%e\t%e\t%e\t%e\t%e\t%e\n", i * dh, laser.x, rc.mask[rc.maskIndex], rc.inputdata[rc.inputIndex], rc.input_signal, laser.x);
						else fprintf(outputfp[0], "%e\t%e\t%e\t%e\t%e\n", i * dh, laser.x, rc.mask[rc.maskIndex], rc.inputdata[rc.inputIndex], rc.input_signal);
						// fprintf(outputfp[0], "%e\t%e\t%e\t%e\t%e\n", i * dh, laser.x, rc.mask[rc.maskIndex], rc.inputdata[rc.inputIndex], rc.input_signal);
					}
					if(i == plotnum) break;
				}
			}
			// InputHistgram(&hist, laser.x);

			if(i % rc.node_interval == rc.node_interval / 2) {
				rc.nodeState[rc.inputIndex - 100][rc.nodeIndex] = laser.x;
				rc.nodeIndex = (rc.nodeIndex + 1) % rc.node_num;
				InputHistgram(&hist, laser.x);
			}
			if(i % rc.mask_interval == 0) rc.maskIndex = (rc.maskIndex + 1) % rc.mask_num;
			if(i % rc.input_interval == 0) rc.inputIndex++;
		}

		end = clock();
		printf("%f\n", (double)(end - start) / CLOCKS_PER_SEC);

		if(out == 1) {
			if(para_flag == 1) fclose(outputfp[0]);
			if(para_flag == 2) {
				// outputfp[0] = fopen("node_state.txt", "w");
				// fprintf(outputfp[0], "Index\tNode\n");
				// for(i = 0; i < rc.allNum; i++) for(j = 0; j < rc.node_num; j++) fprintf(outputfp[0], "%d\t%e\n", i * rc.node_num + j, rc.nodeState[i][j]);
				// fclose(outputfp[0]);

				OutputHistgram(&hist, "histgram.txt\0");
				QR(&rc);
				EvaluateReservoirTrain(&rc);
				outputfp[0] = fopen("rc_training_result.txt", "w");
				fprintf(outputfp[0], "Index\tInput data\tx(t - \\ft\\n\\dD\\n)\tTarget\tTraining result\n");
				for(i = 0; i < rc.trainingNum; i++) fprintf(outputfp[0], "%d\t%e\t%e\t%e\t%e\n", i, rc.inputdata[100 + i], rc.inputdata[100 + i - rc.delay_step] , rc.target[i], rc.output[i]);
				fclose(outputfp[0]);

				EvaluateReservoirTest(&rc);
				outputfp[0] = fopen("rc_test_result.txt", "w");
				fprintf(outputfp[0], "Index\tInput data\tx(t - \\ft\\n\\dD\\n)\tTarget\tTraining result\n");
				for(i = rc.trainingNum; i < rc.allNum; i++) fprintf(outputfp[0], "%d\t%e\t%e\t%e\t%e\n", i, rc.inputdata[100 + i], rc.inputdata[100 + i - rc.delay_step] ,rc.target[i], rc.output[i]);
				fclose(outputfp[0]);
				calc_ser(rc.target, rc.output, rc.testNum);
			}
			break;
		}
		if(out == 2) {
			OutputHistgram(&hist, "histgram.txt\0");
			QR(&rc);
			fprintf(outputfp[0], "%e\t%e\t%e", param, EvaluateReservoirTrain(&rc), EvaluateReservoirTest(&rc));
			fprintf(outputfp[0], "\t%e\t%e\n",calc_ser(rc.target, rc.output, rc.testNum), hist.entropy);

			if(OUT_NODE_STATE == 1) {
				sprintf(filename, ".\\node_state\\node_state_%03d.txt\0", s - amp_min_int);
				outputfp[1] = fopen(filename, "w");
				for(i = 0; i < rc.allNum; i++) {
					for(j = 0; j < rc.node_num; j++) fprintf(outputfp[1], "%e\t", rc.nodeState[i][j]);
					fprintf(outputfp[1], "\n");
				}
				fclose(outputfp[1]);
			}
		}
	}
	//	MAIN calculation ----------------------

	if(out == 2) fclose(outputfp[0]);

	free_rc_memory(&rc);
	free(laser.xdelay);

	return 0;
}

void free_rc_memory(struct Reservoir *rc) {
	int i;

	free(rc->target);
	free(rc->inputdata);
	free(rc->output);
	free(rc->mask);
	free(rc->weight);
	for(i = 0; i < rc->allNum; i++) free(rc->nodeState[i]);
	free(rc->nodeState);

	return;
}

void SetRC_Param(struct Reservoir *rc, double mask_interval, int mask_num, int mask_node_ratio) {
	int i;

	rc->allNum = 100;
	rc->trainingNum = 75;
	rc->testNum = rc->allNum - rc->trainingNum;

	// data points
	rc->target = (double*)malloc(sizeof(double) * 5000);
	rc->inputdata = (double*)malloc(sizeof(double) * 5000);
	rc->output = (double*)malloc(sizeof(double) * 5000);

	srand(1);
	for(i = 0; i < 5000; i++) rc->inputdata[i] = (double)(rand() % 2);
	
//-------------target function ---------------
	for(i = 0; i < 5000; i++){
		if (rc->inputdata[100 + i] == rc->inputdata[100 + (i - rc->delay_step)]){
			rc->target[i] = 0.0;
		}
		else {
			rc->target[i] = 1.0;
		}
	} 
//-------------target function ---------------

	rc->mask_num = mask_num;
	rc->mask_interval = (int)(mask_interval / h + EPS);
	rc->input_interval = rc->mask_interval * rc->mask_num;

	rc->node_num = rc->mask_num * mask_node_ratio;
	rc->node_interval = rc->mask_interval / mask_node_ratio;

	// Make a mask for reservoir computing
	rc->mask = (double*)malloc(sizeof(double) * rc->mask_num);
	MakeMask(rc);

	if(rc->nodeState == NULL) {
		rc->nodeState = (double**)malloc(sizeof(double*) * rc->allNum);
		for(i = 0; i <  rc->allNum; i++) rc->nodeState[i] = (double*)malloc(sizeof(double) * rc->node_num);
	}
	rc->weight = (double*)malloc(sizeof(double) * rc->node_num);

	return;
}

void MakeMask(struct Reservoir *rc) {
	int i, rand_num;

	init_gen_rand(1234);
	for(i = 0; i < rc->mask_num; i++) {
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

double calc_ser(double *target, double *signal, int num)
{
	int count = 0;
	double tmp;
	
	for(int i; i < num; i++)
	{
		if(signal[i] <= 0.5) tmp = 0.0;
		if(signal[i] > 0.5) tmp = 1.0;
		if(tmp != target[i]) count++;
	}
	
	printf("Symbol error rate = %e (Errors : %d / %d)\n", (double)count / (double)num, count, num);
	
	return (double)count / (double)num;
}

void QR(struct Reservoir *rc) {
	int i, j, k;
	double sum, *yx, ridge = 0.001, *R[rc->node_num], *Q[rc->node_num], *QT[rc->node_num];

	for(i = 0; i < rc->node_num; i++) R[i] = (double*)malloc(sizeof(double) * rc->node_num);
	for(i = 0; i < rc->node_num; i++) Q[i] = (double*)malloc(sizeof(double) * rc->node_num);
	for(i = 0; i < rc->node_num; i++) QT[i] = (double*)malloc(sizeof(double) * rc->node_num);
	yx = (double*)malloc(sizeof(double) * rc->node_num);

	for(i = 0; i < rc->node_num; i++) {
		for(j = 0; j < rc->node_num; j++) {
			Q[i][j] = 0.0;
			for(k = 0; k < rc->trainingNum; k++) Q[i][j] += rc->nodeState[k][i] * rc->nodeState[k][j];
			Q[i][j] /= (double)rc->trainingNum;
			R[i][j] = 0.0;
		}
		yx[i] = 0.0;
		for(j = 0; j < rc->trainingNum; j++) yx[i] += rc->target[j] * rc->nodeState[j][i];
		yx[i] /= (double)rc->trainingNum;
	}

	// for(i = 0; i < rc->node_num; i++) Q[i][i] += ridge;/* For ridge regression */

	for(i = 0; i < rc->node_num; i++) {
		for(j = 0; j < i; j++) {
			for(k = 0, R[j][i] = 0.0; k < rc->node_num; k++) R[j][i] += Q[k][i] * Q[k][j];
			for(k = 0; k < rc->node_num; k++) Q[k][i] -= Q[k][j] * R[j][i];
		}
		for(j = 0, R[i][i] = 0.0; j < rc->node_num; j++) R[i][i] += Q[j][i] * Q[j][i];
		R[i][i] = sqrt(R[i][i]);
		for(j = 0; j < rc->node_num; j++) Q[j][i] /= R[i][i];
	}

	for(i = 0; i < rc->node_num; i++) for(j = 0; j < rc->node_num; j++) QT[i][j] = Q[j][i];
	for(i = rc->node_num - 1; i >= 0; i--) {
		for(j = 0, sum = 0.0; j < rc->node_num; j++) sum += QT[i][j] * yx[j];
		for(j = i + 1; j < rc->node_num; j++) sum -= R[i][j] * rc->weight[j];
		rc->weight[i] = sum / R[i][i];
	}

	for(i = 0; i < rc->node_num; i++) free(R[i]);
	for(i = 0; i < rc->node_num; i++) free(Q[i]);
	for(i = 0; i < rc->node_num; i++) free(QT[i]);
	free(yx);

	return;
}

double EvaluateReservoirTrain(struct Reservoir *rc) {
	int i, j;
	double NMSE, sum_y, sqsum_y, variance_y, one_y;

	for(i = 0, NMSE = sum_y = sqsum_y = 0.0; i < rc->trainingNum; i++) {
		sum_y += rc->target[i];
		sqsum_y += rc->target[i] * rc->target[i];
		for(j = 0, one_y = 0.0; j < rc->node_num; j++) one_y += rc->weight[j] * rc->nodeState[i][j];
		NMSE += (rc->target[i] - one_y) * (rc->target[i] - one_y);
		rc->output[i] = one_y;
	}

	variance_y = (sqsum_y - sum_y * sum_y / (double)rc->trainingNum) / (double)rc->trainingNum;
	NMSE /= variance_y * (double)rc->trainingNum;
	printf("Training result: Variance = %e, NMSE = %e\n", variance_y, NMSE);

	return NMSE;
}

double EvaluateReservoirTest(struct Reservoir *rc) {
	int i, j;
	double NMSE, sum_y, sqsum_y, variance_y, one_y;

	for(i = rc->trainingNum, NMSE = sum_y = sqsum_y = 0.0; i < rc->allNum; i++) {
		sum_y += rc->target[i];
		sqsum_y += rc->target[i] * rc->target[i];
		for(j = 0, one_y = 0.0; j < rc->node_num; j++) one_y += rc->weight[j] * rc->nodeState[i][j];
		NMSE += (rc->target[i] - one_y) * (rc->target[i] - one_y);
		rc->output[i] = one_y;
	}

	variance_y = (sqsum_y - sum_y * sum_y / (double)(rc->allNum - rc->trainingNum)) / (double)(rc->allNum - rc->trainingNum);
	NMSE /= variance_y * (double)(rc->allNum - rc->trainingNum);
	printf("Test result: Variance = %e, NMSE = %e\n", variance_y, NMSE);

	return NMSE;
}

void CalcLaserEOfeedback(struct LaserEOfeedback *laser, double input_signal) {
	int i, j, k;
	double a[2], ini[2], b[2][4], cosTheta, feedback, interval_L = h / laser->p.tau_L, interval_H = h / laser->p.tau_H, buffer, beta_mod;

	ini[0] = laser->x;
	ini[1] = laser->y;
	for(i = 0, buffer = 0.0; i < laser->fbNum; i++) buffer += laser->xdelay[laser->fbIndex[i]] * laser->p.weight_fb[i];

	// cosTheta = cos(laser->p.fb_str * buffer + laser->p.phi + input_signal);
	// feedback = laser->p.beta * cosTheta * cosTheta;

	/***** for beta modulation *****/
	cosTheta = cos(laser->p.fb_str * buffer + laser->p.phi);
	beta_mod = 2.0 * cos(input_signal * M_PI / 4.0 - M_PI / 4.0) * cos(input_signal * M_PI / 4.0 - M_PI / 4.0);
	feedback = laser->p.beta * cosTheta * cosTheta * beta_mod;
	/***** for beta modulation *****/

	for(i = 0; i < 4; i++) {
		for(j = 0; j < 2; j++) {
			if(i == 0) a[j] = ini[j];
			else if(i == 1) a[j] = ini[j] + b[j][0] / 2.0;
			else if(i == 2) a[j] = ini[j] + b[j][1] / 2.0;
			else if(i == 3) a[j] = ini[j] + b[j][2];
		}
		b[0][i] = interval_L * ( - (1.0 + laser->p.epsilon) * a[0] - a[1] + feedback);
		b[1][i] = interval_H * a[0];
	}

	for(i = 0; i < 2; i++) a[i] = ini[i] + (b[i][0] + 2.0 * b[i][1] + 2.0 * b[i][2] + b[i][3]) / 6.0;

	laser->xdelay[laser->fbIndex[laser->long_fb]] = laser->x;
	for(i = 0; i < laser->fbNum; i++) laser->fbIndex[i] = (laser->fbIndex[i] + 1) % laser->descretized_num[laser->long_fb];

	a[0] += sqrt( - 2.0 * WGSD * h * log(genrand_open_open()) ) * cos(2.0 * M_PI * genrand_open_open());

	laser->x = a[0];
	laser->y = a[1];

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
	histgram->average = histgram->var = 0.0;
	histgram->n = 0;

	return;
}

void InputHistgram(struct Histdata *histgram, double inputData) {

	if(inputData >= histgram->min && inputData < histgram->max) histgram->data[(int)((inputData - histgram->min) / histgram->div + EPS)]++;
	if(inputData < histgram->minData) histgram->minData = inputData;
	if(inputData > histgram->maxData) histgram->maxData = inputData;
	histgram->average += inputData;
	histgram->var += inputData * inputData;
	histgram->n++;

	return;
}

void OutputHistgram(struct Histdata *histgram, char *filename) {
	int i, maxDatanum;
	FILE *outputfp;

	histgram->average /= histgram->n;
	histgram->var = histgram->var / histgram->n - histgram->average * histgram->average;
	for(i = 0, maxDatanum = 0, histgram->entropy = 0.0; i < histgram->bin; i++) {
		if(histgram->data[i] > maxDatanum) {
			maxDatanum = histgram->data[i];
			histgram->median = (i + 0.5) * histgram->div + histgram->min;
		}
		if(histgram->data[i] != 0.0) histgram->entropy -= (double)histgram->data[i] / (double)(histgram->n) * log2((double)histgram->data[i] / (double)(histgram->n));
	}
	histgram->entropy /= -log2(1.0 / (double)histgram->bin);
	outputfp = fopen(filename, "w");
	fprintf(outputfp, "Maximum data : %e\nMinimum data : %e\nMedian : %e\nAverage : %e\nVariance : %e\nNormalized entropy : %e\n", histgram->maxData, histgram->minData, histgram->median, histgram->average, histgram->var, histgram->entropy);
	fprintf(outputfp, "Data\tProbability\tThe number of points\n");
	for(i = 0; i < histgram->bin; i++) fprintf(outputfp, "%e\t%e\t%d\n", (i + 0.5) * histgram->div + histgram->min, (double)histgram->data[i] / (double)(histgram->n), histgram->data[i]);
	free(histgram->data);
	fclose(outputfp);

	return;
}
