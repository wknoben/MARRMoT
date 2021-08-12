#### m_18_simhyd_7p_3s

comment out lines 83-92 of `m_18_simhyd_7p_3s.m` and substitute with:

```matlab
flux_Ei   = min(S1/delta_t, Ep);
flux_EXC  = P.*(1-(1./(1+exp((S1-insc+.05*insc)/(.01*((.01*insc==0)+insc))))));
flux_INF  = min(coeff.*exp((-1*sq*S2)./smsc),flux_EXC);
flux_INT  = sub*S2/smsc*flux_INF;
flux_REC  = crak*S2/smsc*(flux_INF-flux_INT);
flux_Et   = min([10*S2/smsc,Ep,S2/delta_t]);
flux_SMF  = flux_INF-flux_INT-flux_REC;
flux_GWF  = flux_SMF.*(1-(1./(1+exp((S2-smsc+.05*smsc)/(.01*((.01*smsc==0)+smsc))))));
flux_BAS  = k.*S3;
flux_SRUN = flux_EXC - flux_INF; 
flux_Qt   = flux_SRUN + flux_INT + flux_BAS;
```



#### m_33_sacramento_11p_5s

comment out lines 132-151 of `m_33_sacramento_11p_5s.m` and substitute with:

```matlab
flux_qdir    = pctim.*P;
flux_peff    = (1-pctim).*P;
flux_ru      = ((S2.*uztwm-S1.*uzfwm)/(uztwm+uzfwm))./...
    (1+exp(((S1./uztwm)-(S2./uzfwm)+.05*(S2./uzfwm))/(.01*(((.01*S2./uzfwm) == 0) + (S2./uzfwm)))));
flux_euztw   = min(S1./uztwm.*Ep,S1/delta_t);
flux_twexu   = flux_peff.*(1 - 1./(1+exp((S1-uztwm+.05*uztwm)/(.01*(((.01*uztwm) == 0) + uztwm)))));
flux_qsur    = flux_twexu.*(1 - 1./(1+exp((S2-uzfwm+.05*uzfwm)/(.01*(((.01*uzfwm) == 0) + uzfwm)))));
flux_qint    = kuz.*S2;
flux_euzfw   = min(S2/delta_t,max(0,Ep-flux_euztw));
flux_pc      = max(0,min(S2/delta_t,max(0,S2./uzfwm).*...
    (pbase.*(1+zperc.*((max(0,lztwm-S3)+max(0,lzfwpm-S4)+max(0,lzfwsm-S5))./(lztwm+lzfwpm+lzfwsm)).^(1+rexp)))));
flux_pctw    = (1-pfree) .* flux_pc;
flux_elztw   = min(S3./lztwm.*max(0,Ep-flux_euztw-flux_euzfw),S3/delta_t);
flux_twexl   = flux_pctw.*(1-1./(1+exp((S3-lztwm+.05*lztwm)/(.01*(((.01*lztwm) == 0) + lztwm)))));

if S4 ~= lzfwpm || S5 ~= lzfwsm
    tmp_rd4 = (S4-lzfwpm)/lzfwpm;
    tmp_rd5 = (S5-lzfwsm)/lzfwsm;
    tmp_distr4 = tmp_rd4/(tmp_rd4+tmp_rd5);
    tmp_distr5 = tmp_rd5/(tmp_rd4+tmp_rd5);                
else
    tmp_distr4 = lzfwpm/(lzfwpm+lzfwsm);
    tmp_distr5 = lzfwsm/(lzfwpm+lzfwsm);
end            
flux_twexlp  = tmp_distr4.*flux_twexl;
flux_twexls  = tmp_distr5.*flux_twexl;
flux_pcfwp   = (pfree*tmp_distr4).*flux_pc;
flux_pcfws   = (pfree*tmp_distr5).*flux_pc;
tmp_max = (S4+S5)./(lzfwpm+lzfwsm);
tmp_thresh = 1 ./ (1+exp(((S3./lztwm)-(tmp_max)+.05*(tmp_max))/(.01*(((.01*tmp_max)==0) + tmp_max))));
tmp_sm = (S3.*(lzfwpm+lzfwsm)+lztwm.*(S4+S5))./((lzfwpm+lzfwsm).*(lztwm+lzfwpm+lzfwsm));
flux_rlp     = (S4.*tmp_sm).* tmp_thresh;
flux_rls     = (S5.*tmp_sm).* tmp_thresh;  
flux_qbfp    = klzp.*S4;
flux_qbfs    = klzs.*S5;
```



#### m_37_hbv_15p_5s

comment out lines 117-129 of `m_37_hbv_15p_5s.m` and substitute with:

```matlab
flux_sf   = min(P,max(0,P.*(tt+0.5*tti-T)/tti));
flux_refr = min(max(0,cfr*cfmax*(ttm-T)),S2/delta_t);
flux_melt = min(max(cfmax*(T-ttm),0),S1/delta_t);
flux_rf   = min(P,max(0,P.*(T-(tt-0.5*tti))/tti));
flux_in   = (flux_rf+flux_melt).*(1-...
    (1./(1+exp((S2-whc*S1+.05*whc*S1)/(.01*((.01*whc*S1==0)+(whc*S1)))))));
flux_se   = max((obj.Sold(2)-whc*S1)/delta_t,0);
flux_cf   = min(cflux.*(1-S3/fc),S4/delta_t);
flux_ea   = min([S3/(lp*fc)*Ep,Ep,S3/delta_t]);
flux_r    = (flux_in+flux_se)*((max(S3,0)/fc)^beta);
flux_q0   = min(k0*max(S4,0)^(1+alpha),max(S4/delta_t,0));
flux_perc = min(perc,S4/delta_t);
flux_q1   = k1.*S5;
flux_qt   = uh(1) * (flux_q0 + flux_q1) + stf(1);
```

#### 