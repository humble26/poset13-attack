import Poset13Defs

def checks4 : Bool := (List.range (2^6)).all (fun m => checkAllOn 4 m)
theorem checks4_true : checks4 = true := by native_decide
#print axioms checks4_true
