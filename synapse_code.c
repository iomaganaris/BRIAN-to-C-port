
//an baloume ekswteriki eisodo thelei i<N+1 logika
for (i = 0; i < N; i++){
	if (SpikeArray[i] > 0){
		for (j==0; j < N; j++){
			//auto borei na ginei pio apodotika alla diskoleuei, gi arxi trwme xwro me ta struct synapses
			if (Synapses[i][j].conn){
				Synapses[i][j].FFp = Synapses[i][j].FFP * exp(-(-Synapses[i][j].lastupdate + t)/tau_FFp);
				Synapses[i][j].FBn = Synapses[i][j].FBn * exp(-(-Synapses[i][j].lastupdate + t)/tau_FBn);
				Synapses[i][j].u = Synapses[i][j].U + (-Synapses[i][j].U + Synapses[i][j].u) * exp(-(-Synapses[i][j].lastupdate + t)/tau_u);
				Synapses[i][j].FBp = Synapses[i][j].FBp * exp(-(-Synapses[i][j].lastupdate + t)/tau_FBp);
				Synapses[i][j].R = (Synapses[i][j].R - 1) * exp(-(-Synapse[i][j].lastupdate + t)/tau_r) + 1;
				Synapses[i][j].target_I = s * Synapses[i][j].A * Synapses[i][j].R * Synapses[i][j].u;
				Synapses[i][j].U = Synapses[i][j].U + etaU * (-AFBn * Synapses[i][j].FBn * Synapses[i][j].FBp + AFBp * Synapses[i][j].FBp * Synapses[i][j].FFp);
				if (Synapses[i][j].U < Umin) Synapses[i][j].U = Umin;
				else if (Synapses[i][j].U > Umax) Synapses[i][j].U = Umax;
				Synapses[i][j].w = Synapses[i][j].U * Synapses[i][j].A;
				Synapses[i][j].FFp += 1;
				Synapses[i][j].R -= Synapses[i][j].R * Synapses[i][j].u;
				Synapses[i][j].u += Synapses[i][j].U * (1 - Synapses[i][j].u);
				Synapses[i][j].lastupdate = t;
			}
			//auto theloume? na ekteleitai to ena meta to allo an tixei?
			if (Synapses[j][i].conn){
				Synapses[j][i].FFp = Synapses[j][i].FFp * exp(-(-Synapses[j][i].lastupdate + t)/tau_FFp);
				Synapses[j][i].FBn = Synapses[j][i].FBn * exp(-(-Synapses[j][i].lastupdate + t)/tau_FBn);
				Synapses[j][i].u = Synapses[j][i].U + (-Synapses[j][i].U + Synapses[i][j].u) * exp(-(-Synapses[j][i].lastupdate + t)/tau_u);
				Synapses[j][i].FBp = Synapses[j][i].FBp * exp(-(-Synapses[j][i].lastupdate + t)/tau_FBp);
				Synapses[j][i].R = (Synapses[j][i].R - 1) * exp(-(-Synapse[j][i].lastupdate + t)/tau_r) + 1;
				Synapses[j][i].A = Synapses[j][i].A + etaA * (AFFp * Synapses[j][i].FFp * Synapses[j][i].FBn);
				Synapses[j][i].A = Synapses[j][i].A - etaA * 0.5 * //edw mallon prepei na bgalw meso oro apo kati. to psaxnoume akoma
				if (Synapses[j][i].A < Amin) Synapses[j][i].A = Amin;
				else if (Synapses[j][i].A > Amax) Synapses[j][i].A = Amax;
				Synapses[j][i].w = Synapses[j][i].U * Synapses[j][i].A;
				Synapses[j][i].FBp += 1;
				Synapses[j][i].FBn += 1;
			}
		}
	}
}





typedef struct {
    double conn;
    double w;        // if w = 0 no connection maybe (?)
    double FFp;
    double FBp;
    double FBn;
    double R;
    double u;
    double U;
    double A;
    double target_I;
    double lastupdate;
} Synapse;