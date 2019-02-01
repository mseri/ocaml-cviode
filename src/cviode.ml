(** {2:lib cviode library} *)

type timespan = float * float
(** Representation of a time span. *)

let check _a _f y0 =
  (* TODO: check shapes of a and f, make sure they are all compatible *)
  let _, elts = Owl.Mat.shape y0 in
  assert (Owl.Maths.is_even elts)

let steps t0 t1 dt =
  (t1 -. t0)/.dt |> Float.floor |> int_of_float

(* opening Owl.Mat or Owl.Arr messes up badly all the integer operations *)
let (.${}) = Owl.Mat.(.${})
let (.${}<-) = Owl.Mat.(.${}<-)

let contact1_damped ~a ~f y0 (t0, t1) dt =
  check a f y0;
  let steps = steps t0 t1 dt in
  let _, elts = Owl.Mat.shape y0 in
  let sol = Owl.Mat.empty steps elts in
  sol.${[[0]]}<- y0;
  for idx = 1 to steps-1 do
    let xs = sol.${[[idx-1]; [0; elts/2-1]]} in
    let ps = sol.${[[idx-1]; [elts/2; elts-1]]} in
    let t = t0 +. dt *. (float_of_int idx) in
    let a = a t in
    let c0 = 1.0 -. dt*.a in
    let c1 = 0.5*.dt in
    let fxs = f xs t in
    let xsnew = Owl.Mat.(xs + mul_scalar ps (dt*.c0) + mul_scalar fxs (dt*.c1)) in
    let fxsnew = f xsnew t in
    sol.${[[idx]; [0; elts/2-1]]}<- xsnew;
    sol.${[[idx]; [elts/2; elts-1]]}<- Owl.Mat.(mul_scalar ps c0 + mul_scalar (fxs + fxsnew) c1);
  done;
  sol


let contact2_damped ~a ~f y0 (t0, t1) dt =
  check a f y0;
  let steps = steps t0 t1 dt in
  let _, elts = Owl.Mat.shape y0 in
  let sol = Owl.Mat.empty steps elts in
  sol.${[[0]]}<- y0;
  for idx = 1 to steps-1 do
    let xs = sol.${[[idx-1]; [0; elts/2-1]]} in
    let ps = sol.${[[idx-1]; [elts/2; elts-1]]} in
    let t = t0 +. dt *. (float_of_int idx) in
    let at = a t in
    let fxs = f xs t in
    let c0 = 1.0 -. 0.5*.dt*.at in
    let c0plus = 1.0 +. 0.5*.dt*.at in
    let c1 = 0.5*.dt in
    let xsnew = Owl.Mat.(xs + mul_scalar ps (dt*.c0) + mul_scalar fxs (dt*.c1)) in
    let fxsnew = f xsnew t in
    sol.${[[idx]; [0; elts/2-1]]}<- xsnew;
    sol.${[[idx]; [elts/2; elts-1]]}<- Owl.Mat.(mul_scalar ps (c0/.c0plus) + mul_scalar (fxs + fxsnew) (c1/.c0plus));
  done;
  sol

(* TODO: Add integrator for generic g_2(z) as by description *)