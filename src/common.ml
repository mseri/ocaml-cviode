type timespan = float * float
(** Representation of a time span. *)

let check _a _f y0 =
  (* TODO: check shapes of a and f, make sure they are all compatible *)
  let _, elts = Owl.Mat.shape y0 in
  assert (Owl.Maths.is_even elts)

let steps t0 t1 dt =
  (t1 -. t0)/.dt |> Float.floor |> int_of_float

let integrate ~step y0 (t0, t1) dt =
  (* opening Owl.Mat or Owl.Arr messes up badly all the integer operations *)
  let (.${}) = Owl.Mat.(.${}) in
  let (.${}<-) = Owl.Mat.(.${}<-) in
  let (.%{}) = Owl.Arr.(.%{}) in
  let (.%{}<-) = Owl.Arr.(.%{}<-) in
  (* Allow to get the slices indexes for xs and ys *)
  let slices idx elts =
    [[idx]; [0; elts/2-1]], [[idx]; [elts/2; elts-1]]
  in
  let steps = steps t0 t1 dt in
  let _, elts = Owl.Mat.shape y0 in
  let sol = Owl.Mat.empty steps elts in
  (* These are all fixed step algorithms, saving the time array is
     quite overkill as tspan will just be a deterministic linspace... *)
  let tspan = Owl.Arr.empty [|steps|] in
  sol.${[[0]]}<- y0;
  tspan.%{[|0|]}<-t0;
  for idx = 1 to steps-1 do
    let xi, pi = slices (idx-1) elts in
    let xi', pi' = slices idx elts in
    let xs, ps = sol.${xi}, sol.${pi} in
    let t = tspan.%{[|idx-1|]} in
    let xs', ps', t' = step xs ps t in
    sol.${xi'}<- xs';
    sol.${pi'}<- ps';
    tspan.%{[|idx|]}<-t';
  done;
  tspan, sol