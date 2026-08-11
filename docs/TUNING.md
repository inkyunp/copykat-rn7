# 파라미터 조정 — 어디서 무엇을

**중요: 파라미터 조정은 소스 코드 수정이 아니다.** `ngene.chr`, `LOW.DR`, `UP.DR`,
`KS.cut`, `win.size`는 전부 `copykat()`의 **호출 시 인자**다(포크 소스에는 이 이름들이
함수 시그니처와 변수 참조로만 등장하고, 하드코딩된 상수가 아니다 — `docs/MODIFICATIONS.md` 참조).
따라서 값을 바꿔도 **포크/이미지를 다시 만들 필요가 없다.** `run_copykat_rat.R`에 플래그로
넘기기만 하면 된다.

## 어디서

`run_copykat_rat.R`(이미지 안에서 실행). 예:

```bash
singularity exec --cleanenv --containall --bind "$PWD":/work copykat-rn7.sif \
  Rscript /work/run_copykat_rat.R \
    --input /work/rat_copykat_test_500.rds \
    --output /work/rat_copykat_result \
    --norm-cells /work/norm_cell_names.txt \
    --ngene-chr 5 --low-dr 0.05 --up-dr 0.10 --ks-cut 0.10 --win-size 25
```

같은 이미지에서, 값만 바꿔 여러 번 돌리면 된다.

## 무엇을 (파라미터별 의미와 조정 방향)

| 플래그 | copykat 인자 | 기본 | 무엇을 바꾸나 | 언제 조정 |
| --- | --- | --- | --- | --- |
| `--ngene-chr` | `ngene.chr` | 5 | 각 세포가 **모든 염색체 라벨(rat=21)** 마다 최소 몇 개 유전자를 가져야 통과하는지 | 저depth 세포가 대량 필터링되면 낮춘다(예: 3). rat은 human(23)/mouse(20)과 라벨 수가 달라 필터 강도가 다르다 |
| `--low-dr` | `LOW.DR` | 0.05 | smoothing/segmentation에 쓸 유전자의 **최소 검출 세포 비율** | 데이터가 얕으면 낮춰 유전자를 더 남긴다. `nrow(rawmat)<7000`이면 코드가 `UP.DR<-LOW.DR`로 강제한다 |
| `--up-dr` | `UP.DR` | 0.10 | 최종 결과에 남길 유전자의 최소 검출 비율 | LOW.DR보다 크게. 노이즈가 많으면 올린다 |
| `--ks-cut` | `KS.cut` | 0.10 | segmentation 민감도(KS 검정 컷오프) | **낮추면 더 잘게 분절**(민감), 높이면 보수적. "too few breakpoints"가 뜨면 낮춘다 |
| `--win-size` | `win.size` | 25 | segmentation window당 최소 유전자 수 | 유전자가 적으면 낮춘다. 너무 작으면 노이즈 세그먼트 증가 |

`--cores`(n.cores), `--id-type`(auto/S/E), `--norm-cells`도 조정 가능하다.
`--id-type auto`는 symbol↔`mgi_symbol`, Ensembl↔`ensembl_gene_id` 교집합이 큰 쪽을 고른다.

## 이 데이터(HN00283643, 310셀)에서 관찰된 신호와 조정 힌트

- 로그에 `WARNING: low data quality; assigned LOW.DR to UP.DR`는 뜨지 않았고(8,845 genes past LOW.DR),
  `step 5: segmentation`은 정상 진행됐다. → 필터는 과하지 않았다.
- heatmap이 뚜렷한 arm-level 블록보다 산발적이라면 신호가 약한 것이다. 시도해볼 순서:
  1. `--ks-cut 0.05`로 낮춰 분절 민감도 ↑
  2. 셀 수를 늘려(≥500~1000) 재실행 — 파일럿 310셀은 baseline/clustering이 얕다
  3. `--norm-cells`로 정상 reference를 더 확실히 지정(현재 alveolar macrophage 100셀)
- `ngene.chr`/`LOW.DR`/`UP.DR`는 **실데이터 depth 기준**으로만 정한다. 지금 기본값을 고정 권장값으로
  박지 않는다(데이터마다 다름).

## 소스 수정이 필요한 경우 (참고 — 지금은 아님)

아래는 인자가 아니라 코드/데이터라서, 바꾸려면 포크 수정 + 이미지 재빌드가 필요하다:

- `set.seed(1234)` (copykat 내부 고정 시드)
- 220kb bin 변환을 rat용으로 새로 만드는 것(현재는 mm10처럼 gene-space로 끝냄 — 별도 대공사)
- `full.anno.rn7`의 염색체 포함 범위(현재 1–20+X; Y/MT 포함하려면 `build_full_anno_rn7.R` 수정 후 재빌드)

이번 요청 범위(파라미터 조정)는 전부 인자라서 **이미지 재빌드가 필요 없다.**
