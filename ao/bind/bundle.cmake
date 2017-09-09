file(READ ${CMAKE_CURRENT_LIST_DIR}/ao-guile.cpp AO_GUILE)
file(READ ${CMAKE_CURRENT_LIST_DIR}/shapes.scm SHAPES)
file(READ ${CMAKE_CURRENT_LIST_DIR}/csg.scm CSG)
file(READ ${CMAKE_CURRENT_LIST_DIR}/transforms.scm TRANSFORMS)
file(READ ${CMAKE_CURRENT_LIST_DIR}/vec.scm VEC)

string(REPLACE "AO_GUILE_SHAPES" "${SHAPES}" AO_GUILE "${AO_GUILE}")
string(REPLACE "AO_GUILE_CSG" "${CSG}" AO_GUILE "${AO_GUILE}")
string(REPLACE "AO_GUILE_TRANSFORMS" "${TRANSFORMS}" AO_GUILE "${AO_GUILE}")
string(REPLACE "AO_GUILE_VEC" "${VEC}" AO_GUILE "${AO_GUILE}")

file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/bundle.cpp "${AO_GUILE}")
