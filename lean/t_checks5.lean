import Poset13Defs

def checks5 : Bool := (List.range (2^10)).all (fun m => checkAllOn 5 m)
theorem checks5_true : checks5 = true := by native_decide
#print axioms checks5_true
