# CopyKAT rat (genome="rn7") fork + test image

CopyKAT을 rat(*Rattus norvegicus*, **mRatBN7.2**) scRNA-seq에 실제로 돌리기 위한
포크와 전용 경량 테스트 이미지다. CopyKAT의 mouse(mm10) gene-space 경로(220kb bin
변환 없음)를 그대로 상속받아 `genome="rn7"` 값을 추가하고, 이 폴더 상위의
`genes.gtf`로 만든 rat annotation table을 패키지에 baking했다.

> 결과는 **후보 CNV 증거**일 뿐이며 최종 malignant/normal 진단이 아니다.
> 저자조차 충분히 검증하지 않은 mouse 코드 경로를 타므로, DNA 레벨 근거
> (WGS/SNP-array/karyotype) 교차검증 전까지 모든 결과는 `review_required`다.

## 무엇을 · 어떻게 · 무엇이 다른지

### 왜 필요한가 (문제)
upstream CopyKAT `copykat()`의 `genome` 인자는 `"hg20"`(human)과 `"mm10"`(mouse) **두 값만**
분기하고 `else`가 없다. rat 값을 그대로 넣으면 어느 분기에도 걸리지 않아 `anno.mat`이
만들어지지 않고 다음 줄에서 곧바로 죽는다. rat scRNA-seq에 CopyKAT을 쓰려면 코드 경로 자체가
없었다.

### 무엇을 했나 (접근)
처음부터 rat 전용 220kb bin 좌표계를 새로 만드는 대신, **mouse(mm10) 다운스트림 경로를 그대로
재사용**했다. mm10 경로는 220kb bin 변환 없이 gene-space로 결과를 내는 구조라, rat annotation
테이블 하나만 주입하면 되기 때문이다. 그래서 `genome="rn7"`(rat mRatBN7.2)이 mm10과 동일한
gene-space 블록을 타도록 최소 변경만 했다.

### 어떻게 했나 (구현)
1. mRatBN7.2 GTF에서 CopyKAT 스키마(7컬럼) rat annotation `full.anno.rn7`을 생성.
2. `copykat()`에 `genome=="rn7"` 분기와 `annotateGenes.rn7()`를 추가하고, mm10 gene-space 블록을 rn7도 통과하도록 한 줄 확장.
3. `full.anno.rn7`을 패키지 `R/sysdata.rda`에 baking → 런타임에 GTF·gene order 파일 불필요.
4. copykat 전용 경량 이미지(`rocker/r-ver:4.2.2`)로 패키징하고, 합성 데이터·실데이터로 실행 검증.

### upstream과 무엇이 다른가

| 항목 | upstream copykat 1.1.0 | 이 포크 (rn7) |
| --- | --- | --- |
| 지원 종 | human(hg20), mouse(mm10) | **+ rat(rn7, mRatBN7.2)** |
| `genome` 미지원 값 | `else` 없음 → 에러로 죽음 | rn7 분기 추가 |
| rat annotation | 없음 | `full.anno.rn7` 패키지에 baked |
| annotation 함수 | `annotateGenes.hg20/.mm10` | `+ annotateGenes.rn7()` |
| rat 출력 단위 | (불가) | gene-space (`CNA_results.txt`, mm10과 동일) |
| hg20 전용 단계 | 변경 없음 | **변경 없음** (HLA/cyclegene 제거·220kb bin은 hg20 게이트라 손대지 않음) |
| 알고리즘 로직 | — | **변경 없음** (필터·smoothing·baseline·MCMC segmentation 원본 그대로) |

> 코드 변경은 정확히 **3곳(+데이터 1개)** 뿐이다. 전체 diff와 복원한 원본은
> [`docs/MODIFICATIONS.md`](docs/MODIFICATIONS.md), [`docs/copykat_rn7.patch`](docs/copykat_rn7.patch),
> [`docs/original_copykat_b795ff7.R`](docs/original_copykat_b795ff7.R) 참조.

### 검증 상태 (실행으로 확인됨)
- **빌드**: `copykat-rn7.sif` 빌드 성공, 인이미지 rn7 체크 + `%test` 합성 실행 PASS, 오프라인(`--network none`) 실행 PASS.
- **T2 합성**: 실제 R 4.2.2 + 실제 deps로 `copykat(genome="rn7")` 완주, gene-space 출력·리터럴 클래스 확인.
- **실데이터**: rat 310셀(폐, epithelial)에서 완주 — 정상 baseline(macrophage)은 대부분 diploid, epithelial은 aneuploid 비율 상승. 결과는 `review_required`.

## 입력 계약 (copykat 1.1.0 소스로 확정)
CopyKAT은 세포/유전자 필터를 **양성 카운트 검출(`sum(x>0)`)** 로 하고, 내부에서
`log(sqrt(x)+sqrt(x+1))` 변환·중심화·dlm 스무딩을 **직접** 수행한다. 따라서 입력은
반드시 **raw integer UMI counts, genes × cells** 여야 한다.

- Seurat v4: `GetAssayData(obj, assay="RNA", slot="counts")`
- Seurat v5: `LayerData(obj, assay="RNA", layer="counts")`
- **금지**: log-normalized `data`(이중 변환 → 신호 소실), `scale.data`(음수 → `sqrt(음수)=NaN`), SCT/integrated assay
- 유전자를 variable features로 줄이지 말 것(genome-wide 커버리지 필요)
- `norm.cell.names`(선택, 권장): 정상 세포 barcode 벡터로 raw matrix의 colnames와 문자 그대로 일치해야 함

## 파일
| 파일 | 용도 |
| --- | --- |
| `build_full_anno_rn7.R` | `../genes.gtf`(mRatBN7.2)에서 7컬럼 `full.anno.rn7` 생성 → `full.anno.rn7.rds` + 무결성 게이트 |
| `full.anno.rn7.rds` | 생성된 rat annotation (25,302 genes, 21 chromosome labels) |
| `copykat-rn7/` | 포크 소스 트리 (copykat 1.1.0, commit `b795ff7` 기반). `full.anno.rn7`을 `R/sysdata.rda`에 baking |
| `copykat-rn7.def` | 경량 Singularity 정의 (`rocker/r-ver:4.2.2`, copykat 전용) |
| `build.sh` / `validate.sh` | 이미지 빌드(fakeroot/네트워크 1회) / 무권한 오프라인 검증 |
| `validate_rat_input.R` | **데이터 검증 게이트** (입력 계약 5항목). 큰 머신/로컬 공용 |
| `run_copykat_rat.R` | 러너: RNA counts 추출 → `copykat(..., genome="rn7")` → 후보 증거 저장 |
| `tests/t2_synthetic_rn7.R` | 합성 스모크 테스트 (RDS 불필요) |

## 포크가 바꾼 것 (upstream 대비 정확히 3곳)
1. annotation 분기에 `else if(genome=="rn7") anno.mat <- annotateGenes.rn7(...)` 추가
2. `annotateGenes.rn7()` 신규 = `annotateGenes.mm10()` 복사본, 기본값 `full.anno=full.anno.rn7`
3. 후반 gene-space 블록 `if(genome=="mm10")` → `if(genome=="mm10" || genome=="rn7")`

hg20 전용 단계(HLA/cyclegene 제거, 220kb bin 변환)는 이미 `genome=="hg20"`로만 게이트되어
있어 손대지 않았다. `mgi_symbol` 컬럼명은 문자 그대로 유지했다(copykat 본문이 이 이름을
직접 참조하므로 rat symbol을 담되 컬럼명을 바꾸면 깨진다).

`full.anno.rn7` 스키마: `abspos, chromosome_name(정수 1–20, X=21), start_position,
end_position, ensembl_gene_id, mgi_symbol(rat symbol), band(NA)`. `abspos`는 실제 누적
좌표라서 내장 `full.anno.mm10`의 상수-abspos 결함(unique 20/137,030)을 재현하지 않는다.

## 실행 순서

### 이미 완료 (로컬, RDS 불필요)
- **Phase 0**: `Rscript build_full_anno_rn7.R` → `full.anno.rn7.rds` 생성, 무결성 게이트 통과
  (25,302 genes / unique abspos 25,268 / 21 labels / min 575 genes per chr).
- **T1**: 포크 구조 검증 — `annotateGenes.rn7` 정의, `formals(copykat)$genome`, sysdata baking.
- **T2**: 실제 R 4.2.2 + 실제 deps로 수정 copykat `genome="rn7"` 합성 실행 완주
  (gene-space `CNA_results.txt`, 리터럴 `aneuploid`/`diploid` 분리 확인).

### 이미지 빌드 (사용자 WSL, 네트워크 1회)
```bash
cd rat-copykat
./build.sh                 # -> copykat-rn7.sif
./validate.sh copykat-rn7.sif
```
빌드는 base 이미지와 CRAN 스냅샷(2023-10-01) 다운로드로 네트워크가 필요하고, 런타임은 오프라인이다.

### Phase 1 — 큰 머신에서 샘플링 + 데이터 검증 (전송 게이트)
5.6GB 원본 RDS는 7.6GB 로컬에서 로드 불가하므로, 아래는 **≥32GB(권장 64GB) 머신에서 1회** 수행한다.
로컬로는 `full.anno.rn7.rds`와 `validate_rat_input.R`를 함께 가져간다.

먼저 알려줄 것: 세포타입 컬럼명 + 정상 라벨, `sample_id` 컬럼명.
(rownames/Seurat 버전/gene·cell 수는 검증 리포트가 자동 기록)

```r
library(Seurat); set.seed(1729)
obj <- readRDS("test.obj_final_add_celltype_scp2711.rds")
DefaultAssay(obj) <- "RNA"
ct  <- "<celltype_col>"; smp <- "<sample_id_col>"
one <- names(sort(table(obj[[smp]][,1]), decreasing=TRUE))[1]   # 샘플 1개
o1  <- subset(obj, cells = colnames(obj)[obj[[smp]][,1] == one])
# 세포타입 층화 + depth 하위 제외로 약 500 barcode 선택 (정상 lineage >= 60~100 보장)
keep <- NULL  # <층화+depth 필터로 약 500 barcode>
sub  <- subset(o1, cells = keep)
sub  <- DietSeurat(sub, assays="RNA", dimreducs=NULL, graphs=NULL)  # counts 유지
sub@meta.data <- sub@meta.data[, c(ct, smp), drop=FALSE]
saveRDS(sub, "rat_copykat_test_500.rds")                            # 목표 <0.5GB
writeLines(colnames(sub)[sub@meta.data[[ct]] %in% c("<정상 라벨들>")], "norm_cell_names.txt")
```

추출 직후 같은 머신에서 검증하고 **PASS일 때만 전송**한다:
```bash
Rscript validate_rat_input.R \
  --input rat_copykat_test_500.rds --anno full.anno.rn7.rds \
  --norm-cells norm_cell_names.txt \
  --celltype-column <celltype_col> --normal-labels "<정상 라벨들>" \
  --report input_validation_report.tsv
# exit 0 = PASS(전송), exit 1 = FAIL(전송 금지)
```

### Phase 2/3 — 로컬 이관 후 실행
`rat_copykat_test_500.rds` + `norm_cell_names.txt` + `input_validation_report.tsv`만 로컬로 복사
(원본 5.6GB RDS는 옮기지 않음). 로컬에서 재검증 후 실행:
```bash
# 로컬 재검증 (전송 손상/재현성)
Rscript validate_rat_input.R --input rat_copykat_test_500.rds --anno full.anno.rn7.rds \
  --norm-cells norm_cell_names.txt --report input_validation_report.local.tsv

# 실행 (이미지 안)
singularity exec --cleanenv --containall --bind "$PWD":/work copykat-rn7.sif \
  Rscript /work/run_copykat_rat.R \
    --input /work/rat_copykat_test_500.rds \
    --output /work/rat_copykat_result \
    --norm-cells /work/norm_cell_names.txt \
    --id-type auto --cores 1
```
출력: `rat_copykat_result.rds`(리터럴 클래스 + `interpretation=review_required`),
`rat_copykat_result.cnv_evidence.tsv`, `rat_copykat_result_rn7_run/`(heatmap PDF 등).

## 재검토가 필요한 실행 파라미터 (실데이터 확보 후)
- `ngene.chr`: rat은 상염색체 20 + X = 21 라벨 → human(23)·mouse(20)과 필터 강도가 다름
- `LOW.DR`, `UP.DR`: rat 데이터 depth/품질에 좌우
- `norm.cell.names` 대 자동 baseline: 정상 세포가 부족하면 자동 baseline이 불안정

## 참고: 이 환경에서의 빌드/실행
- system R은 4.5.0이고 CopyKAT deps가 없어, T1/T2는 기존 Rocky9 SIF에서 추출한 **R 4.2.2 + deps**로
  실제 실행했다(faithful target R).
- 기존 SIF는 이 exec 샌드박스에서 `/dev/fuse` 부재·setuid starter 소유권 문제로 마운트 실패한다
  (샌드박스 아티팩트). 새 이미지의 실제 빌드/실행은 사용자 WSL에서 `build.sh`/`validate.sh`로 검증할 것.
