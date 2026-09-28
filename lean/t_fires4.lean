import Poset13Defs

theorem cFires4 : cFires 4 (List.range (2^6)) > 0 := by native_decide
theorem pFires4 : pFires 4 (List.range (2^6)) > 0 := by native_decide
#print axioms cFires4
#print axioms pFires4
