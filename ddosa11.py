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
        ht['GODOL']=self.input_ecorrdata.godol
        ht['mcecDOL']=self.input_ecorrdata.mcecdol
        ht['l2reDOL']=self.input_ecorrdata.l2redol
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

