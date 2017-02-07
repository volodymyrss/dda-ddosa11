from ddosa import *

class ibis_isgr_energy(DataAnalysis):
    cached=False

    input_scw=ScWData
    input_ecorrdata=GetEcorrCalDB

    version="v6_extras"

    def main(self):

        remove_withtemplate("isgri_events_corrected.fits(ISGR-EVTS-COR.tpl)")

        construct_gnrl_scwg_grp(self.input_scw,[\
            self.input_scw.scwpath+"/isgri_events.fits[3]", \
            self.input_scw.scwpath+"/ibis_hk.fits[IBIS-DPE.-CNV]", \
            self.input_scw.auxadppath+"/time_correlation.fits[AUXL-TCOR-HIS]" 
        ])

        import_attr(self.input_scw.scwpath+"/swg.fits",['OBTSTART','OBTEND'])
        set_attr({'ISDCLEVL':"COR"})

        bin="ibis_isgr_energy"
        ht=heatool(bin)
        ht['inGRP']="og.fits"
        ht['outCorEvts']="isgri_events_corrected.fits(ISGR-EVTS-COR.tpl)"
        ht['useGTI']="y"
        #ht['eraseALL']="y"
        ht['randSeed']=500
        ht['riseDOL']="auto"
        #ht['riseDOL']=self.input_ecorrdata.risedol
        ht['GODOL']="auto"#self.input_ecorrdata.godol
        ht['mcecDOL']="auto"#self.input_ecorrdata.mcecdol
        ht['l2reDOL']="auto"#self.input_ecorrdata.l2redol
        ht['chatter']="10"
        ht.run()


        self.output_events=DataFile("isgri_events_corrected.fits")

class ibis_comp_energy(DataAnalysis):
    cached=False

    input_scw=ScWData
    input_ecorrdata=GetEcorrCalDB

    version="v6_extras"

    def main(self):

        remove_withtemplate("compton_events_corrected.fits(COMP-SGLE-COR.tpl)")

        construct_gnrl_scwg_grp(self.input_scw,[\
            self.input_scw.scwpath+"/compton_events.fits[3]", \
            self.input_scw.scwpath+"/ibis_hk.fits[IBIS-DPE.-CNV]", \
            self.input_scw.auxadppath+"/time_correlation.fits[AUXL-TCOR-HIS]" 
        ])

        import_attr(self.input_scw.scwpath+"/swg.fits",['OBTSTART','OBTEND'])
        set_attr({'ISDCLEVL':"COR"})

        bin="ibis_comp_energy"
        ht=heatool(bin)
        ht['inGRP']="og.fits"
        ht['outGRP']=""
        ht['outCorEvts']="compton_events_corrected.fits(COMP-SGLE-COR.tpl)"
        ht['useGTI']="y"
        ht['extName']="COMP-SGLE-RAW"
        #ht['eraseALL']="y"
        ht['randSeed']=500
       # ht['riseDOL']="auto"
        ht['riseDOL']="auto"
        ht['GODOL']=self.input_ecorrdata.godol
        ht['mcecDOL']=self.input_ecorrdata.mcecdol
        ht['l2reDOL']=self.input_ecorrdata.l2redol
        ht['enerDOL']="auto"
        ht['chatter']=5
        ht.run()


        self.output_events=DataFile("compton_events_corrected.fits")

#/sps/integral/data/ic/ic_snapshot_20140321/ic/ibis/mod/isgr_effc_mod_0001.fits

class GetEcorrCalDB(DataAnalysis):
    input=["ecorr_standard_OSA10.2"]
    input_lut2=GetLUT2
    input_ibisic=IBIS_ICRoot
    input_scw=ScWData

    cached=False

    #ignore_input=["input"]


    def main(self):
        self.godol=self.input_ibisic.ibisicroot+"/cal/ibis_isgr_gain_offset_0010.fits"
        self.risedol=self.input_ibisic.ibisicroot+"/mod/isgr_rise_mod_0239.fits"
        self.mcecdol=self.input_ibisic.ibisicroot+"/mod/isgr_mcec_mod_0001.fits"
        self.l2redol=self.input_ibisic.ibisicroot+"/mod/isgr_l2re_mod_0001.fits"


class BinEventsVirtual(DataAnalysis):
    input_scw=ScWData

    input_ibisic=IBIS_ICRoot

    input_events=ISGRIEvents
    input_gti=ibis_gti
    input_dead=ibis_dead

    target_level=None
    input_bins=None

    maxrisetime=116
    minrisetime=16

    version="v3"
    
    cached=True

    default_log_level="binevents"

    ii_shadow_build_binary="ii_shadow_build"


    def get_version(self):
        v=self.get_signature()+"."+self.version
        if self.maxrisetime!=116:
            v+=".lrt%i"%self.maxrisetime
        if self.minrisetime!=16:
            v+=".hrt%i"%self.minrisetime
        return v

    def main(self):
        if self.target_level is None or self.input_bins is None:
            raise Exception("VirtualAnalysis: please inherit!")
        
        self.pre_process()

        # ask stephane why need raw events
        construct_gnrl_scwg_grp(self.input_scw,[\
            self.input_events.events.get_path(), \
            self.input_scw.scwpath+"/isgri_events.fits[3]", \
            self.input_scw.scwpath+"/ibis_hk.fits[IBIS-DPE.-CNV]", \
            self.input_scw.auxadppath+"/time_correlation.fits[AUXL-TCOR-HIS]" \
        ]) # get separately tc etc

        import_attr(self.input_scw.scwpath+"/swg.fits",['OBTSTART','OBTEND'])
        set_attr({'ISDCLEVL':self.target_level})

        det_fn="isgri_detector_shadowgram_%s.fits"%self.target_level
        det_tpl="(ISGR-DETE-SHD-IDX.tpl)"
        eff_fn="isgri_efficiency_shadowgram_%s.fits"%self.target_level
        eff_tpl="(ISGR-EFFI-SHD-IDX.tpl)"

        remove_withtemplate(det_fn+det_tpl)
        remove_withtemplate(eff_fn+eff_tpl)

        bin=self.ii_shadow_build_binary
        ht=heatool(bin)
        ht['outSWGGRP']="og.fits[GROUPING,1,BINTABLE]"
        ht['inDead']=self.input_dead.output_dead.get_path()
        ht['inGTI']=self.input_gti.output_gti.get_path()
        ht['gti_name'] = 'MERGED_ISGRI'
        ht['outputLevel'] = self.target_level

        print "target_level",self.target_level

        print "has rmfbins "+str(self.input_bins.rmfbins) if hasattr(self.input_bins,'rmfbins') else "no rmfbins"
        print "has binrmfext" if hasattr(self.input_bins,'binrmfext') else "no binrmfext"

        if ( self.target_level=="BIN_I" or not hasattr(self.input_bins,'rmfbins') or not self.input_bins.rmfbins or not hasattr(self.input_bins,'binrmfext') ) and not ( hasattr(self.input_bins,'rmfbins') and self.input_bins.rmfbins ): # fix!!
            ht['isgri_e_num'] = len(self.input_bins.bins)
            ht['isgri_e_min'] = " ".join([str(a[0]) for a in self.input_bins.bins])
            ht['isgri_e_max'] = " ".join([str(a[1]) for a in self.input_bins.bins])
        elif self.target_level=="BIN_S" or ( hasattr(self.input_bins,'rmfbins') and self.input_bins.rmfbins ) or hasattr(self.input_bins,'binrmfext'):
            ht['isgri_e_num'] = -1

            rmf=self.input_bins.binrmfext
            if isinstance(rmf,DataFile):
                ht['inEnergyValues'] = rmf.get_path()  #  +"[1]" will this work?
            else:
                ht['inEnergyValues'] = rmf #  +"[1]" will this work?
        else:
            raise Exception("neither bins given!")

        ht['isgri_min_rise'] = self.minrisetime
        ht['isgri_max_rise'] = self.maxrisetime
        ht['isgri_t_len'] = 10000000
        ht['idxLowThre']=self.input_scw.revdirpath+"/idx/isgri_context_index.fits[1]"
        ht['idxNoisy']=self.input_scw.revdirpath+"/idx/isgri_prp_noise_index.fits[1]"
        ht['outRawShadow']=det_fn+det_tpl
        ht['outEffShadow']=eff_fn+eff_tpl
        ht['inEFFC']="auto"
        #ht['inEFFC']=self.input_ibisic.ibisicroot+"/mod/isgr_effc_mod_0001.fits"

        self.extra_pars(ht)

        ht.run()

        self.shadow_detector=DataFile(det_fn)
        self.shadow_efficiency=DataFile(eff_fn)
        
        self.post_process()

    def extra_pars(selt,ht):
        pass

    def post_process(self):
        pass
    
    def pre_process(self):
        pass



class BinEventsImage(BinEventsVirtual):
    target_level="BIN_I"
    input_bins=ImageBins

class BinEventsSpectra(BinEventsVirtual):
    target_level="BIN_S"
    input_bins=SpectraBins

class SpectraBins(DataAnalysis):
    input_binsname="spectral_bins_256"
    
    rmfbins=True

    version="v1"
    def main(self):
        self.binrmf=os.environ['CURRENT_IC']+"/ic/ibis/mod/isgr_ebds_mod_0001.fits"
        e=pyfits.open(self.binrmf)[1].data
        self.bins=zip(e['E_MIN'],e['E_MAX'])
        self.binrmfext=self.binrmf+'[1]'

    def get_binrmfext(self):
        return self.binrmfext


