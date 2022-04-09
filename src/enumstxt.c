/**
 * \file enumstxt.c
 * \brief Source file to define text strings for enumerated data types.
 * \author see AUTHORS.
 * \copyright see AUTHORS.
 */

#include <glib.h>

#include "text.h"
#include "enumstxt.h"

const char *NodeTxt[] = { t_JUNCTION,
  t_RESERVOIR,
  t_TANK
};

const char *LinkTxt[] = { w_CV,
  w_PIPE,
  w_PUMP,
  w_PRV,
  w_PSV,
  w_PBV,
  w_FCV,
  w_TCV,
  w_GPV
};

const char *StatTxt[] = { t_XHEAD,
  t_TEMPCLOSED,
  t_CLOSED,
  t_OPEN,
  t_ACTIVE,
  t_XFLOW,
  t_XFCV,
  t_XPRESSURE,
  t_FILLING,
  t_EMPTYING,
  t_OVERFLOWING
};

const char *FormTxt[] = { w_HW,
  w_DW,
  w_CM
};

const char *RptFormTxt[] = { t_HW,
  t_DW,
  t_CM
};

const char *RptFlowUnitsTxt[] = { u_CFS,
  u_GPM,
  u_MGD,
  u_IMGD,
  u_AFD,
  u_LPS,
  u_LPM,
  u_MLD,
  u_CMH,
  u_CMD
};

const char *FlowUnitsTxt[] = { w_CFS,
  w_GPM,
  w_MGD,
  w_IMGD,
  w_AFD,
  w_LPS,
  w_LPM,
  w_MLD,
  w_CMH,
  w_CMD
};

const char *PressUnitsTxt[] = { w_PSI,
  w_KPA,
  w_METERS
};

const char *DemandModelTxt[] = { w_DDA,
  w_PDA,
  NULL
};

const char *SourceTxt[] = { w_CONCEN,
  w_MASS,
  w_SETPOINT,
  w_FLOWPACED
};

const char *ControlTxt[] = { w_BELOW,
  w_ABOVE,
  w_TIME,
  w_CLOCKTIME
};

const char *TstatTxt[] = { w_NONE,
  w_AVG,
  w_MIN,
  w_MAX,
  w_RANGE
};

const char *MixTxt[] = { w_MIXED,
  w_2COMP,
  w_FIFO,
  w_LIFO,
  NULL
};

const char *RptFlagTxt[] = { w_NO,
  w_YES,
  w_FULL
};

const char *SectTxt[] = { s_TITLE, s_JUNCTIONS, s_RESERVOIRS,
  s_TANKS, s_PIPES, s_PUMPS,
  s_VALVES, s_CONTROLS, s_RULES,
  s_DEMANDS, s_SOURCES, s_EMITTERS,
  s_PATTERNS, s_CURVES, s_QUALITY,
  s_STATUS, s_ROUGHNESS, s_ENERGY,
  s_REACTIONS, s_MIXING, s_REPORT,
  s_TIMES, s_OPTIONS, s_COORDS,
  s_VERTICES, s_LABELS, s_BACKDROP,
  s_TAGS, s_END,
  NULL
};

const char *Fldname[] = { t_ELEV, t_DEMAND, t_HEAD,
  t_PRESSURE, t_QUALITY, t_LENGTH,
  t_DIAM, t_FLOW, t_VELOCITY,
  t_HEADLOSS, t_LINKQUAL, t_LINKSTATUS,
  t_SETTING, t_REACTRATE, t_FRICTION,
  "", "", "", "", "", "", NULL
};
