/* externals for DOOM */

#ifdef _LINUX_
#ifdef _PGF77_
#define doom_blow       doom_blow_
#define doom_close      doom_close_
#define doom_compress   doom_compress_
#define doom_comp_max   doom_comp_max_
#define doom_debug      doom_debug_
#define doom_dump       doom_dump_
#define doom_dyncode    doom_dyncode_
#define doom_exist      doom_exist_
#define doom_exist_env  doom_exist_env_
#define doom_fetch      doom_fetch_
#define doom_galign     doom_galign_
#define doom_gbase      doom_gbase_
#define doom_gcorrect   doom_gcorrect_
#define doom_gdirect    doom_gdirect_
#define doom_gelcode    doom_gelcode_
#define doom_gelement   doom_gelement_
#define doom_getenv     doom_getenv_
#define doom_getvar     doom_getvar_
#define doom_gfelem     doom_gfelem_
#define doom_gfield     doom_gfield_
#define doom_gline      doom_gline_
#define doom_gmonitor   doom_gmonitor_
#define doom_gthead     doom_gthead_
#define doom_gtbody     doom_gtbody_
#define doom_gparam     doom_gparam_
#define doom_gparent    doom_gparent_
#define doom_gpos       doom_gpos_
#define doom_gsequ      doom_gsequ_
#define doom_gstring    doom_gstring_
#define doom_gsumm      doom_gsumm_
#define doom_gsurv      doom_gsurv_
#define doom_gtime      doom_gtime_
#define doom_gtime_env  doom_gtime_env_
#define doom_gtrack     doom_gtrack_
#define doom_gturn      doom_gturn_
#define doom_gtype      doom_gtype_
#define doom_list       doom_list_
#define doom_lpos       doom_lpos_
#define doom_memory     doom_memory_
#define doom_nadd       doom_nadd_
#define doom_newdb      doom_newdb_
#define doom_npos       doom_npos_
#define doom_open       doom_open_
#define doom_palign     doom_palign_
#define doom_pcorrect   doom_pcorrect_
#define doom_pdirect    doom_pdirect_
#define doom_pelement   doom_pelement_
#define doom_pfelem     doom_pfelem_
#define doom_pfield     doom_pfield_
#define doom_pline      doom_pline_
#define doom_plink      doom_plink_
#define doom_pmonitor   doom_pmonitor_
#define doom_pthead     doom_pthead_
#define doom_gtnames    doom_gtnames_
#define doom_ptbody     doom_ptbody_
#define doom_pparam     doom_pparam_
#define doom_pparent    doom_pparent_
#define doom_prtab      doom_prtab_
#define doom_psequ      doom_psequ_
#define doom_pstring    doom_pstring_
#define doom_psumm      doom_psumm_
#define doom_psurv      doom_psurv_
#define doom_ptrack     doom_ptrack_
#define doom_pturn      doom_pturn_
#define doom_ptype      doom_ptype_
#define doom_recall     doom_recall_
#define doom_rline      doom_rline_
#define doom_robj       doom_robj_
#define doom_save       doom_save_
#define doom_setenv     doom_setenv_
#define doom_setvar     doom_setvar_
#define doom_sobj       doom_sobj_
#define doom_snoup      doom_snoup_
#define doom_store      doom_store_
#define doom_swap       doom_swap_
#define doom_time       doom_time_
#define doom_wtree      doom_wtree_
#endif
#endif
#ifdef _LINUX_
#ifndef _PGF77_
#define doom_blow       doom_blow__
#define doom_close      doom_close__
#define doom_compress   doom_compress__
#define doom_comp_max   doom_comp_max__
#define doom_debug      doom_debug__
#define doom_dump       doom_dump__
#define doom_dyncode    doom_dyncode__
#define doom_exist      doom_exist__
#define doom_exist_env  doom_exist_env__
#define doom_fetch      doom_fetch__
#define doom_galign     doom_galign__
#define doom_gbase      doom_gbase__
#define doom_gcorrect   doom_gcorrect__
#define doom_gdirect    doom_gdirect__
#define doom_gelcode    doom_gelcode__
#define doom_gelement   doom_gelement__
#define doom_getenv     doom_getenv__
#define doom_getvar     doom_getvar__
#define doom_gfelem     doom_gfelem__
#define doom_gfield     doom_gfield__
#define doom_gline      doom_gline__
#define doom_gmonitor   doom_gmonitor__
#define doom_gthead     doom_gthead__
#define doom_gtnames    doom_gtnames__
#define doom_gtbody     doom_gtbody__
#define doom_gparam     doom_gparam__
#define doom_gparent    doom_gparent__
#define doom_gpos       doom_gpos__
#define doom_gsequ      doom_gsequ__
#define doom_gstring    doom_gstring__
#define doom_gsumm      doom_gsumm__
#define doom_gsurv      doom_gsurv__
#define doom_gtime      doom_gtime__
#define doom_gtime_env  doom_gtime_env__
#define doom_gtrack     doom_gtrack__
#define doom_gturn      doom_gturn__
#define doom_gtype      doom_gtype__
#define doom_list       doom_list__
#define doom_lpos       doom_lpos__
#define doom_memory     doom_memory__
#define doom_nadd       doom_nadd__
#define doom_newdb      doom_newdb__
#define doom_npos       doom_npos__
#define doom_open       doom_open__
#define doom_palign     doom_palign__
#define doom_pcorrect   doom_pcorrect__
#define doom_pdirect    doom_pdirect__
#define doom_pelement   doom_pelement__
#define doom_pfelem     doom_pfelem__
#define doom_pfield     doom_pfield__
#define doom_pline      doom_pline__
#define doom_plink      doom_plink__
#define doom_pmonitor   doom_pmonitor__
#define doom_pthead     doom_pthead__
#define doom_ptbody     doom_ptbody__
#define doom_pparam     doom_pparam__
#define doom_pparent    doom_pparent__
#define doom_prtab      doom_prtab__
#define doom_psequ      doom_psequ__
#define doom_pstring    doom_pstring__
#define doom_psumm      doom_psumm__
#define doom_psurv      doom_psurv__
#define doom_ptrack     doom_ptrack__
#define doom_pturn      doom_pturn__
#define doom_ptype      doom_ptype__
#define doom_recall     doom_recall__
#define doom_rline      doom_rline__
#define doom_robj       doom_robj__
#define doom_save       doom_save__
#define doom_setenv     doom_setenv__
#define doom_setvar     doom_setvar__
#define doom_sobj       doom_sobj__
#define doom_snoup      doom_snoup__
#define doom_store      doom_store__
#define doom_swap       doom_swap__
#define doom_time       doom_time__
#define doom_wtree      doom_wtree__
#endif
#endif
#ifndef _LINUX_
#define doom_blow       doom_blow_
#define doom_close      doom_close_
#define doom_compress   doom_compress_
#define doom_comp_max   doom_comp_max_
#define doom_debug      doom_debug_
#define doom_dump       doom_dump_
#define doom_dyncode    doom_dyncode_
#define doom_exist      doom_exist_
#define doom_exist_env  doom_exist_env_
#define doom_fetch      doom_fetch_
#define doom_galign     doom_galign_
#define doom_gbase      doom_gbase_
#define doom_gcorrect   doom_gcorrect_
#define doom_gdirect    doom_gdirect_
#define doom_gelcode    doom_gelcode_
#define doom_gelement   doom_gelement_
#define doom_getenv     doom_getenv_
#define doom_getvar     doom_getvar_
#define doom_gfelem     doom_gfelem_
#define doom_gfield     doom_gfield_
#define doom_gline      doom_gline_
#define doom_gmonitor   doom_gmonitor_
#define doom_gthead     doom_gthead_
#define doom_gtbody     doom_gtbody_
#define doom_gparam     doom_gparam_
#define doom_gparent    doom_gparent_
#define doom_gpos       doom_gpos_
#define doom_gsequ      doom_gsequ_
#define doom_gstring    doom_gstring_
#define doom_gsumm      doom_gsumm_
#define doom_gsurv      doom_gsurv_
#define doom_gtime      doom_gtime_
#define doom_gtime_env  doom_gtime_env_
#define doom_gtrack     doom_gtrack_
#define doom_gturn      doom_gturn_
#define doom_gtype      doom_gtype_
#define doom_list       doom_list_
#define doom_lpos       doom_lpos_
#define doom_memory     doom_memory_
#define doom_nadd       doom_nadd_
#define doom_newdb      doom_newdb_
#define doom_npos       doom_npos_
#define doom_open       doom_open_
#define doom_palign     doom_palign_
#define doom_pcorrect   doom_pcorrect_
#define doom_pdirect    doom_pdirect_
#define doom_pelement   doom_pelement_
#define doom_pfelem     doom_pfelem_
#define doom_pfield     doom_pfield_
#define doom_pline      doom_pline_
#define doom_plink      doom_plink_
#define doom_pmonitor   doom_pmonitor_
#define doom_pthead     doom_pthead_
#define doom_gtnames    doom_gtnames_
#define doom_ptbody     doom_ptbody_
#define doom_pparam     doom_pparam_
#define doom_pparent    doom_pparent_
#define doom_prtab      doom_prtab_
#define doom_psequ      doom_psequ_
#define doom_pstring    doom_pstring_
#define doom_psumm      doom_psumm_
#define doom_psurv      doom_psurv_
#define doom_ptrack     doom_ptrack_
#define doom_pturn      doom_pturn_
#define doom_ptype      doom_ptype_
#define doom_recall     doom_recall_
#define doom_rline      doom_rline_
#define doom_robj       doom_robj_
#define doom_save       doom_save_
#define doom_setenv     doom_setenv_
#define doom_setvar     doom_setvar_
#define doom_sobj       doom_sobj_
#define doom_snoup      doom_snoup_
#define doom_store      doom_store_
#define doom_swap       doom_swap_
#define doom_time       doom_time_
#define doom_wtree      doom_wtree_
#endif

/* define DOOM special table codes */

#define ALIGN_CODE         101
#define FIELD_CODE         102
#define TRACK_SV           103
#define TRACK_RECORD       104
#define CORR_CODE          105
#define SEQU_CODE          106
#define POLY_CODE          107
#define ELPAR_CODE         108
#define TABLE_HEAD         109
#define TABLE_BODY         110
#define TABLE_CODE         111
#define ENTRY_CODE         112
#define EXIT_CODE          113
#define APERTURE_CODE      114
#define MON_CODE           115
#define MAP_CODE           116

/* define overall parameters */

#define CONT_INT 10
#define CONT_OBJ 28
#define C_DUM_SIZE 4000
#define DOUBLE_SIZE 8
#define DYNAMIC_CODE 200
#define INFO_BASE 20
#define INFO_STEP 10
#define INT_SIZE 4
#define KEY_LENGTH 48
#define LINE_LENGTH 80
#define N_SYSTEM 5
#define O_LIST_SIZE 40000
#define REC_CONT 4
#define TABLE_OFFSET 2
#define MAX_DELTA 25
#define OPTICS_SUMM 17
#define OPTICS_INTL 4

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

  /* variables */

extern const int el_att_num;
extern const int el_att_pos[];
extern char* el_att_names[];

/* Fortran callable routines  */

extern void doom_blow(char*, int*, int*, int*);
extern void doom_close();
extern void doom_compress(char*, int*, int*);
extern void doom_comp_max(char*, int*, int*, int*);
extern int  doom_debug();
extern void doom_dump(int*);
extern int  doom_dyncode();
extern int  doom_exist(char*, int*);
extern int  doom_exist_env(char*, char*, int*);
extern void db_get_info();
extern void db_get_system();
extern void db_update();
extern struct object* doom_fetch(char*);
extern void doom_galign(char*, int*, int*, double*);
extern void doom_gbase(char*, int*, int*, char*, int*);
extern void doom_gcorrect(char*, int*, int*, double*);
extern void doom_gdirect(char*, int*, int*, char*);
extern int  doom_gelcode(char*);
extern void doom_gelement(char*, int*, double*);
extern void doom_getenv(char*, char*);
extern void doom_gfelem(char*, int*, int*, int*, double*, char*);
extern void doom_gfield(char*, int*, int*, double*);
extern void doom_gline(char*, char*, int*);
extern void doom_gmonitor(char*, int*, int*, double*);
extern void doom_gthead(char*, int*, char*, int*, int*, int*, int*, double*);
extern void doom_gtnames(char*, int*, int*, char*);
extern void doom_gtbody(char*, double*, char*, int*, int*, int*, char*, 
                        int*, double*);
extern void doom_gparam(char*, int*, double*, int*, char*);
extern void doom_gparent(char*, int*, int*, char*, int*);
extern int  doom_gpos(char*);
extern void doom_gsequ(char*, int*, int*, char*, int*, double*);
extern void doom_gstring(char*, int*, int*, char*);
extern void doom_gsumm(char*, int*, double*);
extern void doom_gsurv(char*, int*, int*, char*, int*, double*);
extern void doom_gtime(char*, int*, double*);
extern void doom_gtime_env(char*, char*, int*, double*);
extern void doom_gtrack(int*, int*);
extern void doom_gturn(int*, int*, int*, double*);
extern void doom_gtype(char*, int*, int*, char*, int*);
extern void doom_list();
extern int  doom_lpos(int*, int*, int*);
extern void doom_memory();
extern int  doom_nadd(char*, char*, int*, int*);
extern int  doom_newdb();
extern int  doom_npos(char*, char*, int*, int*);
extern void doom_open(char*);
extern void doom_palign(char*, int*, int*, double*);
extern void doom_pcorrect(char*, int*, int*, double*);
extern void doom_pelement(char*, int*, double*);
extern void doom_pfelem(char*, int*, char*, int*, int*, int*, double*, char*);
extern void doom_pdirect(char*, int*, int*, char*);
extern void doom_pfield(char*, int*, int*, double*);
extern void doom_pline(char*, char*, int*);
extern void doom_plink(char*, int*, int*);
extern void doom_pmonitor(char*, int*, int*, double*);
extern void doom_pthead(char*, char*, int*, int*, int*, double*,
                         int*, int*, char*, int*);
extern void doom_ptbody(char*, char*, int*, int*, int*, 
                        char*, int*, double*);
extern void doom_pparam(char*, int*, double*, int*, char*);
extern void doom_pparent(char*, int*, char*);
extern void doom_prtab(char*, int*, int*, char*, int*);
extern void doom_psequ(char*, int*, int*, char*, int*, double*);
extern void doom_pstring(char*, int*, int*, char*);
extern void doom_psumm(char*, int*, double*);
extern void doom_psurv(char*, int*, int*, char*, int*, double*);
extern void doom_ptrack(int*, int*);
extern void doom_pturn(int*, int*, int*, double*);
extern void doom_ptype(char*, int*, char*);
extern void doom_recall(char*, int*, int*, int*, int*, 
                int*, double*, char*);
extern void doom_rline(char*);
extern void doom_robj(char*, int*, char*, char*,
                      int*, int*, int*, int*, int*, double*, char*);
extern void doom_save(struct object*);
extern void doom_setenv(char*, char*);
extern void doom_setvar(char*, int*);
extern void doom_sobj(char*, int*, char*, char*,
                      int*, int*, int*, int*, int*, double*, char*);
extern void doom_snoup(int*);
extern void doom_store(char*, int*, int*, int*, int*, 
                int*, double*, char*);
extern void doom_swap(int*);
extern void doom_time(double*);
extern void doom_wtree(int*, int*, int*);

/* internal routines */

extern void clean_obj(struct object*);
extern int db_close();
extern int db_exists(char*);
extern void db_get(char*, int*, int*, int*, double*, char*);
extern void db_get_info();
extern void db_get_system();
extern int db_is_open();
extern int db_length(int*);
extern void db_list();
extern int db_open(char*, int);
extern int db_put(char*, int*, int*, int*, double*, char*);
extern void db_put_info();
extern void db_put_system();
extern void delete_obj(struct object*);
extern void dmp_i_list(int, int*);
extern void error_dump(int*);
extern void fill_obj(struct object*, int, int, int, const int*, 
                     const double*, const char*);
extern void free_obj(struct object*);
extern struct object* get_obj(char*);
extern struct object* get_obj_list();
extern int get_text(int, int, char*);
extern void grow_obj(struct object*, int, int, int, int);
extern int i_list_add(int, struct object*);
extern int i_list_pos(int, struct object*);
extern void make_dollar_key(const char*, const char*, char*);
extern void make_key(const char*, const int*, char*);
extern void make_names(struct object*);
extern struct object* make_obj(char*, int, int, int, int);
extern struct object* make_isp_table();
extern void make_volatile();
extern void mycpy(char*, char*);
extern void myscpy(char*, char*);
extern void my_time();
extern void n_list_add(int, char*, struct object*);
extern int n_list_pos(char*, struct object*);
extern void n_list_replace(int, char*, struct object*);
extern int non_zero(int*, int);
extern void o_list_add(int, struct object*);
extern void open_new();
extern void open_old();
extern void prt_obj(struct object*);
extern void prt_obj_full(struct object*);
extern void prt_i_list(struct object*);
extern void prt_o_list(struct object*);
extern void put_obj(struct object*);
extern int put_text(int, int, char*);
extern void rec_get(int*, int*, int*, double*, char*);
extern void rec_put(int*, int*, int*, double*, char*);
extern void shift_char(char*, int, int);
extern int start_get(char*, int, int*, int*, int*, double*, char*);
extern void swap_int(int, int*);
extern void swap_double(int, int*);
extern void upper_conv(char*);
extern int var_list_pos(char*, struct object*);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif
