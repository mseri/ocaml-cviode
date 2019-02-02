(** {2:lib cviode library} *)

type timespan = float * float
(** Representation of a time span. *)

(** Private function, allow to get the slices indexes for xs and ys *)
let slices idx elts =
  [[idx]; [0; elts/2-1]], [[idx]; [elts/2; elts-1]]

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
    let xi, pi = slices (idx-1) elts in
    let xi', pi' = slices idx elts in
    let xs, ps = sol.${xi}, sol.${pi} in
    let t = t0 +. dt*.(float_of_int idx) in
    let c0 = 1.0 -. dt*.(a t) in
    let c1 = 0.5 *. dt in
    let fxs = f xs t in
    let xs' = Owl.Mat.(xs + ps *$ (dt*.c0) + fxs *$ (dt*.c1)) in
    let fxs' = f xs' t in
    sol.${xi'}<- xs';
    sol.${pi'}<- Owl.Mat.(ps *$ c0 + (fxs + fxs') *$ c1);
  done;
  sol


let contact2_damped ~a ~f y0 (t0, t1) dt =
  check a f y0;
  let steps = steps t0 t1 dt in
  let _, elts = Owl.Mat.shape y0 in
  let sol = Owl.Mat.empty steps elts in
  sol.${[[0]]}<- y0;
  for idx = 1 to steps-1 do
    let xi, pi = slices (idx-1) elts in
    let xi', pi' = slices idx elts in
    let xs, ps = sol.${xi}, sol.${pi} in
    let t = t0 +. dt *. (float_of_int idx) in
    let at = a t in
    let fxs = f xs t in
    let c0m = 1.0 -. 0.5*.dt*.at in
    let c0p = 1.0 +. 0.5*.dt*.at in
    let c1 = 0.5*.dt in
    let xs' = Owl.Mat.(xs + ps *$ (dt*.c0m) + fxs *$ (dt*.c1)) in
    let fxs' = f xs' t in
    sol.${xi'}<- xs';
    sol.${pi'}<- Owl.Mat.(ps *$ (c0m/.c0p) + (fxs + fxs') *$ (c1/.c0p));
  done;
  sol

(* TODO: Add integrator for generic g_2(z) as by description *)