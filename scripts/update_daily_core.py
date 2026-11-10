
      - name: Validate original wind fields
        run: |
          python - <<'PY'
          import json
          import math
          from pathlib import Path

          folder = Path("data/latest")

          with (folder / "manifest.json").open() as f:
              manifest = json.load(f)

          expected = {
              "wind850": 850,
              "wind250": 250,
          }

          assert set(manifest["wind_layers"]) == set(expected)

          for name, level in expected.items():

              with (folder / f"{name}.json").open() as f:
                  data = json.load(f)

              assert data["variable"] == name
              assert data["level_hpa"] == level
              assert data["date"] == manifest["date"]
              assert data["field_type"] == "original_daily_mean_wind"
              assert data["units"] == "m/s"
              assert data["shape"] == [256, 512]

              assert len(data["lat"]) == 256
              assert len(data["lon"]) == 512
              assert len(data["u"]) == 256
              assert len(data["v"]) == 256

              valid = 0

              for urow, vrow in zip(data["u"], data["v"]):

                  assert len(urow) == 512
                  assert len(vrow) == 512

                  for u, v in zip(urow, vrow):

                      assert (u is None) == (v is None)

                      if u is not None:
                          assert math.isfinite(u)
                          assert math.isfinite(v)
                          valid += 1

              assert valid > 0

              print(
                  f"SUCCESS {name}: "
                  f"{valid:,} original wind vectors"
              )

          print("SUCCESS: Both wind layers validated")
          PY
