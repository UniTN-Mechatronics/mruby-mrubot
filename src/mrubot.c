/***************************************************************************/
/*                                                                         */
/* mrubot.c - mruby testing                                                */
/* Copyright (C) 2015 Paolo Bosetti                                        */
/* paolo[dot]bosetti[at]unitn.it                                           */
/* Department of Industrial Engineering, University of Trento              */
/*                                                                         */
/* This library is free software.  You can redistribute it and/or          */
/* modify it under the terms of the GNU GENERAL PUBLIC LICENSE 2.0.        */
/*                                                                         */
/* This library is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* Artistic License 2.0 for more details.                                  */
/*                                                                         */
/* See the file LICENSE                                                    */
/*                                                                         */
/***************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <sys/param.h>
#include <time.h>

#include "mruby.h"
#include "mruby/variable.h"
#include "mruby/string.h"
#include "mruby/data.h"
#include "mruby/class.h"
#include "mruby/value.h"
#include "mruby/array.h"
#include "mruby/numeric.h"

#include "memory.h"

// Struct holding data:
typedef struct {
  double d;
  int i;
  double *ary;
} mrubot_data_s;

// Garbage collector handler, for mrubot_data struct
// if mrubot_data contains other dynamic data, free it too!
// Check it with GC.start
static void mrubot_data_destructor(mrb_state *mrb, void *p_) {
  mrubot_data_s *pd = (mrubot_data_s *)p_;
  free(pd->ary);
  free(pd);
  // or simply:
  // mrb_free(mrb, pd);
};

// Creating data type and reference for GC, in a const struct
const struct mrb_data_type mrubot_data_type = {"mrubot_data",
                                               mrubot_data_destructor};

// Utility function for getting the struct out of the wrapping IV @data
static void mrb_mrubot_get_data(mrb_state *mrb, mrb_value self,
                                mrubot_data_s **data) {
  mrb_value data_value;
  data_value = mrb_iv_get(mrb, self, mrb_intern_lit(mrb, "@data"));

  // Loading data from data_value into p_data:
  Data_Get_Struct(mrb, data_value, &mrubot_data_type, *data);
  if (!*data)
    mrb_raise(mrb, E_RUNTIME_ERROR, "Could not access @data");
}

// Data Initializer C function (not exposed!)
static void mrb_mrubot_init(mrb_state *mrb, mrb_value self, double d) {
  mrb_value data_value;  // this IV holds the data
  mrubot_data_s *p_data; // pointer to the C struct

  data_value = mrb_iv_get(mrb, self, mrb_intern_lit(mrb, "@data"));

  // if @data already exists, free its content:
  if (!mrb_nil_p(data_value)) {
    Data_Get_Struct(mrb, data_value, &mrubot_data_type, p_data);
    free(p_data);
  }
  // Allocate and zero-out the data struct:
  p_data = malloc(sizeof(mrubot_data_s));
  memset(p_data, 0, sizeof(mrubot_data_s));
  if (!p_data)
    mrb_raise(mrb, E_RUNTIME_ERROR, "Could not allocate @data");

  // Wrap struct into @data:
  mrb_iv_set(
      mrb, self, mrb_intern_lit(mrb, "@data"), // set @data
      mrb_obj_value(                           // with value hold in struct
          Data_Wrap_Struct(mrb, mrb->object_class, &mrubot_data_type, p_data)));

  // Now set values into struct:
  p_data->d = d;
  p_data->i = 10;
  p_data->ary = malloc(sizeof(double) * p_data->i);
  memset(p_data->ary, 0, sizeof(double) * p_data->i);
}

static mrb_value mrb_mrubot_initialize(mrb_state *mrb, mrb_value self) {
  mrb_value d = mrb_nil_value();
  mrb_get_args(mrb, "o", &d);

  // Call strcut initializer:
  mrb_mrubot_init(mrb, self, mrb_to_flo(mrb, d));
  return mrb_nil_value();
}

static mrb_value mrb_mrubot_d(mrb_state *mrb, mrb_value self) {
  mrubot_data_s *p_data = NULL;

  // call utility for unwrapping @data into p_data:
  mrb_mrubot_get_data(mrb, self, &p_data);

  // Play with p_data content:
  return mrb_float_value(mrb, p_data->d);
}

static mrb_value mrb_mrubot_set_d(mrb_state *mrb, mrb_value self) {
  mrb_value d_value = mrb_nil_value();
  mrubot_data_s *p_data = NULL;

  mrb_get_args(mrb, "o", &d_value);

  // call utility for unwrapping @data into p_data:
  mrb_mrubot_get_data(mrb, self, &p_data);

  p_data->d = mrb_to_flo(mrb, d_value);
  return d_value;
}

static mrb_value mrb_mrubot_ary(mrb_state *mrb, mrb_value self) {
  mrubot_data_s *p_data = NULL;
  mrb_value ary;
  mrb_int i;
  // call utility for unwrapping @data into p_data:
  mrb_mrubot_get_data(mrb, self, &p_data);

  // Play with p_data content:
  ary = mrb_ary_new_capa(mrb, p_data->i);
  for (i = 0; i < p_data->i; i++) {
    mrb_ary_set(mrb, ary, i, mrb_float_value(mrb, p_data->ary[i]));
  }
  return ary;
}

static mrb_value mrb_mrubot_set_ary(mrb_state *mrb, mrb_value self) {
  mrb_value ary_in = mrb_nil_value();
  mrubot_data_s *p_data = NULL;
  mrb_int i;
  mrb_value elem;
  mrb_get_args(mrb, "A", &ary_in);

  // call utility for unwrapping @data into p_data:
  mrb_mrubot_get_data(mrb, self, &p_data);
  if (p_data->ary)
    free(p_data->ary);

  p_data->i = RARRAY_LEN(ary_in);
  p_data->ary = malloc(sizeof(double) * p_data->i);
  for (i = 0; i < p_data->i; i++) {
    elem = mrb_ary_entry(ary_in, i);
    if (mrb_fixnum_p(elem) || mrb_float_p(elem))
      p_data->ary[i] = mrb_to_flo(mrb, elem);
    else {
      p_data->i = 0;
      p_data->ary = malloc(0);
      mrb_raisef(mrb, E_RUNTIME_ERROR, "Non-numeric entry at position %S",
                 mrb_fixnum_value(i));
    }
  }
  return mrb_fixnum_value(i);
}

/* BALLISTIC */
#define G 9.81

static mrb_value mrb_ballistic_theta_for_xy(mrb_state *mrb, mrb_value self) {
  mrb_float x, y, v, v2, sr, theta1, theta2;
  mrb_value result = mrb_ary_new_capa(mrb, 2);

  x = mrb_to_flo(mrb, mrb_iv_get(mrb, self, mrb_intern_lit(mrb, "@x")));
  y = mrb_to_flo(mrb, mrb_iv_get(mrb, self, mrb_intern_lit(mrb, "@y")));
  v = mrb_to_flo(mrb, mrb_iv_get(mrb, self, mrb_intern_lit(mrb, "@v")));
  v2 = pow(v, 2);
  sr = sqrt(pow(v2, 2) - G * (G * pow(x, 2) + 2 * y * v2));
  theta1 = atan2((v2 - sr), (G * x));
  theta2 = atan2((v2 + sr), (G * x));
  mrb_ary_push(mrb, result, mrb_float_value(mrb, theta1));
  mrb_ary_push(mrb, result, mrb_float_value(mrb, theta2));
  return result;
}

static mrb_value mrb_ballistic_theta_for_x(mrb_state *mrb, mrb_value self) {
  mrb_int nargs;
  mrb_float x;
  nargs = mrb_get_args(mrb, "|f", &x);
  if (nargs == 1)
    mrb_iv_set(mrb, self, mrb_intern_lit(mrb, "@x"), mrb_float_value(mrb, x));
  return mrb_ballistic_theta_for_xy(mrb, self);
}

static mrb_value mrb_ballistic_theta_for_y(mrb_state *mrb, mrb_value self) {
  mrb_int nargs;
  mrb_float y;
  nargs = mrb_get_args(mrb, "|f", &y);
  if (nargs == 1)
    mrb_iv_set(mrb, self, mrb_intern_lit(mrb, "@y"), mrb_float_value(mrb, y));
  return mrb_ballistic_theta_for_xy(mrb, self);
}

static mrb_value mrb_ballistic_reach_angle(mrb_state *mrb, mrb_value self) {
  mrb_float d, v;
  mrb_value result;
  mrb_int nargs;
  nargs = mrb_get_args(mrb, "|f", &d);
  if (nargs < 1)
    d = mrb_to_flo(mrb, mrb_iv_get(mrb, self, mrb_intern_lit(mrb, "@x")));
  v = mrb_to_flo(mrb, mrb_iv_get(mrb, self, mrb_intern_lit(mrb, "@v")));
  result = mrb_float_value(mrb, asin((G * d) / pow(v, 2)) / 2.0);
  return result;
}

static mrb_value mrb_ballistic_distance(mrb_state *mrb, mrb_value self) {
  mrb_float theta, v, dy, result;
  mrb_int nargs;
  nargs = mrb_get_args(mrb, "|ff", &theta, &dy);
  if (nargs == 0) { // assume angle as @theta and dy as 0
    theta =
        mrb_to_flo(mrb, mrb_iv_get(mrb, self, mrb_intern_lit(mrb, "@theta")));
    dy = 0.0;
  }
  if (nargs == 1) { // read angle and assume dy = 0
    dy = 0.0;
  }
  v = mrb_to_flo(mrb, mrb_iv_get(mrb, self, mrb_intern_lit(mrb, "@v")));
  result = (v * cos(theta) / G) *
           (v * sin(theta) + sqrt(pow(v * sin(theta), 2) + 2 * G * dy));
  return mrb_float_value(mrb, result);
}

/* Newton-Raphson */
static mrb_value mrb_newton_solve(mrb_state *mrb, mrb_value self) {
  mrb_float x0, x1, h, y, yd, max_error;
  mrb_int itr, max_iter;
  mrb_value f, fd;
  mrb_value block = mrb_nil_value();
  mrb_value block_args;

  mrb_int nargs = mrb_get_args(mrb, "|f&", &x0, &block);

  if (nargs == 0) {
    x0 = mrb_to_flo(mrb, mrb_iv_get(mrb, self, mrb_intern_lit(mrb, "@x0")));
  }
  max_error =
      mrb_to_flo(mrb, mrb_iv_get(mrb, self, mrb_intern_lit(mrb, "@max_error")));
  max_iter =
      mrb_fixnum(mrb_iv_get(mrb, self, mrb_intern_lit(mrb, "@max_iter")));
  f = mrb_iv_get(mrb, self, mrb_intern_lit(mrb, "@f"));
  fd = mrb_iv_get(mrb, self, mrb_intern_lit(mrb, "@fd"));

  if (!mrb_nil_p(block)) {
    block_args = mrb_ary_new_capa(mrb, 2);
  }
  x1 = 0;
  for (itr = 0; itr < max_iter; itr++) {
    y = mrb_to_flo(mrb,
                   mrb_funcall(mrb, f, "call", 1, mrb_float_value(mrb, x0)));
    yd = mrb_to_flo(mrb,
                    mrb_funcall(mrb, fd, "call", 1, mrb_float_value(mrb, x0)));
    h = y / yd;
    x1 = x0 - h;

    if (!mrb_nil_p(block)) {
      mrb_ary_set(mrb, block_args, 0, mrb_fixnum_value(itr));
      mrb_ary_set(mrb, block_args, 1, mrb_float_value(mrb, x1));
      mrb_ary_set(mrb, block_args, 2, mrb_float_value(mrb, y));
      mrb_yield(mrb, block, block_args);
    }

    if (fabs(h) < max_error) {
      break;
    }
    x0 = x1;
  }
  return mrb_float_value(mrb, x1);
}

/* MEMORY INFO */
static mrb_value mrb_process_getCurrentRSS(mrb_state *mrb, mrb_value self) {
  return mrb_fixnum_value(getCurrentRSS());
}

static mrb_value mrb_process_getPeakRSS(mrb_state *mrb, mrb_value self) {
  return mrb_fixnum_value(getPeakRSS());
}

static mrb_value mrb_kernel_daemon(mrb_state *mrb, mrb_value self) {
  mrb_bool nochdir, noclose;
  mrb_int nargs = mrb_get_args(mrb, "|bb", &nochdir, &noclose);
  if (nargs == 1) {
    noclose = 1;
  } else if (nargs == 0) {
    nochdir = 1;
    noclose = 1;
  }
  if (0 != daemon(nochdir, noclose))
    mrb_raisef(mrb, E_RUNTIME_ERROR, "Could not daemonize.\n%S",
               mrb_str_new_cstr(mrb, strerror(errno)));
  return mrb_true_value();
}

static mrb_value mrb_kernel_daemon(mrb_state *mrb, mrb_value self) {
  char buf[MAXPATHLEN];
  mrb_value result;
  int daemonized;
  mrb_bool nochdir = 0, noclose = 0;
  mrb_get_args(mrb, "|bb", &nochdir, &noclose);
  result = mrb_ary_new_capa(mrb, 2);
  daemonized = daemon(nochdir, noclose);
  getcwd(buf, MAXPATHLEN);
  mrb_ary_push(mrb, result, mrb_fixnum_value(daemonized));
  mrb_ary_push(mrb, result, mrb_str_new_cstr(mrb, buf));
  mrb_gv_set(mrb, mrb_intern_lit(mrb, "$PID"), mrb_fixnum_value(getpid()));
  mrb_gv_set(mrb, mrb_intern_lit(mrb, "$PPID"), mrb_fixnum_value(getppid()));
  return result;
}

static mrb_value mrb_kernel_sleep(mrb_state *mrb, mrb_value self) {
  mrb_float period;
  struct timespec ts = {}, rts = {};
  mrb_get_args(mrb, "f", &period);

  ts.tv_sec = (mrb_int)period;
  ts.tv_nsec = (mrb_int)((period - ts.tv_sec) * 1e9);
  if (0 != nanosleep(&ts, &rts)) {
    double actual = rts.tv_sec + rts.tv_nsec / (double)1e9;
    mrb_value actual_v = mrb_float_value(mrb, actual);
    char *buf = NULL;
    asprintf(&buf, "Sleep interrupted (errno: '%s'). Slept for %f s",
             strerror(errno), actual);
    mrb_value exc =
        mrb_exc_new(mrb, mrb_class_get(mrb, "SleepError"), buf, strlen(buf));
    mrb_iv_set(mrb, exc, mrb_intern_lit(mrb, "@actual"), actual_v);
    free(buf);
    mrb_exc_raise(mrb, exc);
  }
  return mrb_float_value(mrb, 0);
}

void mrb_mruby_mrubot_gem_init(mrb_state *mrb) {
  struct RClass *process_mod, *mrubot_mod;    // Modules
  struct RClass *mrubot, *ballistic, *newton; // Classes
  // Kernel methods:
  mrb_define_method(mrb, mrb->kernel_module, "daemon", mrb_kernel_daemon,
                    MRB_ARGS_OPT(2));
  mrb_define_method(mrb, mrb->kernel_module, "sleep", mrb_kernel_sleep,
                    MRB_ARGS_REQ(1));
  mrb_load_string(mrb,
                  "class SleepError < Exception; attr_reader :actual; end");
  mrb_define_method(mrb, mrb->kernel_module, "daemon", mrb_kernel_daemon,
                    MRB_ARGS_OPT(2));

  mrubot_mod = mrb_define_module(mrb, "Mrubot");

  ballistic =
      mrb_define_class_under(mrb, mrubot_mod, "Ballistic", mrb->object_class);
  mrb_const_set(mrb, mrb_obj_value(ballistic), mrb_intern_lit(mrb, "G"),
                mrb_float_value(mrb, G));
  mrb_define_method(mrb, ballistic, "theta_for_xy", mrb_ballistic_theta_for_xy,
                    MRB_ARGS_NONE());
  mrb_define_method(mrb, ballistic, "theta_for_x", mrb_ballistic_theta_for_x,
                    MRB_ARGS_OPT(1));
  mrb_define_method(mrb, ballistic, "theta_for_y", mrb_ballistic_theta_for_y,
                    MRB_ARGS_OPT(1));
  mrb_define_method(mrb, ballistic, "reach_angle", mrb_ballistic_reach_angle,
                    MRB_ARGS_OPT(1));
  mrb_define_method(mrb, ballistic, "distance", mrb_ballistic_distance,
                    MRB_ARGS_OPT(1));

  newton = mrb_define_class_under(mrb, mrubot_mod, "NewtonRaphson",
                                  mrb->object_class);
  mrb_define_method(mrb, newton, "solve", mrb_newton_solve, MRB_ARGS_NONE());

  mrubot = mrb_define_class_under(mrb, mrubot_mod, "Mrubot", mrb->object_class);
  mrb_define_method(mrb, mrubot, "initialize", mrb_mrubot_initialize,
                    MRB_ARGS_NONE());
  mrb_define_method(mrb, mrubot, "d", mrb_mrubot_d, MRB_ARGS_NONE());
  mrb_define_method(mrb, mrubot, "d=", mrb_mrubot_set_d, MRB_ARGS_REQ(1));
  mrb_define_method(mrb, mrubot, "ary", mrb_mrubot_ary, MRB_ARGS_NONE());
  mrb_define_method(mrb, mrubot, "ary=", mrb_mrubot_set_ary, MRB_ARGS_REQ(1));

  process_mod = mrb_define_module(mrb, "ProcessInfo");
  mrb_const_set(mrb, mrb_obj_value(process_mod), mrb_intern_lit(mrb, "PID"),
                mrb_fixnum_value(getpid()));
  mrb_const_set(mrb, mrb_obj_value(process_mod), mrb_intern_lit(mrb, "PPID"),
                mrb_fixnum_value(getppid()));

  mrb_define_class_method(mrb, process_mod, "current_mem",
                          mrb_process_getCurrentRSS, MRB_ARGS_NONE());
  mrb_define_class_method(mrb, process_mod, "peak_mem", mrb_process_getPeakRSS,
                          MRB_ARGS_NONE());
}

void mrb_mruby_mrubot_gem_final(mrb_state *mrb) {}
