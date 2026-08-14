# Источники Day 4

- `SAND2013-2569`, Table 4 — длина, станции, хорды, twist и x-offset.
- `SAND2013-2569`, Table 5 — Saertex DB и foam: толщины и упругие свойства.
- `SAND2013-2569`, Table 11 — ширина spar-cap region 600 мм.
- `SAND2013-2569`, Table 13 — две DB-обшивки на сторону и foam core 50 мм.
- NREL archived FAST `CertTest/5MW_Baseline/Airfoils` — координаты восьми профилей.
- `Wind_wing/results/structural/glass_station_results.csv` — factored parked-fault
  moment envelope, использованный для нагрузки `V = |dM/dr|`.

`cad_inputs.zip` содержит одну таблицу 19 CAD-станций и восемь официальных
файлов координат профилей. Архив читается напрямую `code/freecad_builder.py`.

Для новых расчётов приняты первичные значения Sandia, а не старая локальная БД:
Saertex DB `t=1 мм, Ex=13,6 ГПа, Ey=13,3 ГПа, Gxy=11,8 ГПа`; foam
`Ex=Ey=256 МПа, Gxy=22 МПа, ρ=200 кг/м³`. Foam shear strength `1,8 МПа`
остаётся явно помеченным консервативным H100 proxy, потому что Table 5 прочность
foam не приводит.
