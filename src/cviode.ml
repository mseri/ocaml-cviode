open Common

(** {2:lib cviode library} *)
let contact1_damped ~a ~f y0 tspan dt =
  check a f y0;
  let step xs ps t0 =
    let t = t0 +. dt in
    let c0 = 1.0 -. dt*.(a t) in
    let c1 = 0.5 *. dt in
    let fxs = f xs t in
    let xs' = Owl.Mat.(xs + ps *$ (dt*.c0) + fxs *$ (dt*.c1)) in
    let fxs' = f xs' t in
    let ps' = Owl.Mat.(ps *$ c0 + (fxs + fxs') *$ c1) in
    xs', ps', t
  in
  integrate ~step y0 tspan dt


let contact2_damped ~a ~f y0 tspan dt =
  check a f y0;
  let step xs ps t0 =
    let t = t0 +. dt in
    let at = a t in
    let fxs = f xs t in
    let c0m = 1.0 -. 0.5*.dt*.at in
    let c0p = 1.0 +. 0.5*.dt*.at in
    let c1 = 0.5*.dt in
    let xs' = Owl.Mat.(xs + ps *$ (dt*.c0m) + fxs *$ (dt*.c1)) in
    let fxs' = f xs' t in
    let ps' = Owl.Mat.(ps *$ (c0m/.c0p) + (fxs + fxs') *$ (c1/.c0p)) in
    xs', ps', t
  in
  integrate ~step y0 tspan dt

(* TODO: Add integrator for generic g_2(z) as by description *)