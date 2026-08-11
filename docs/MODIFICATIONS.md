# CopyKAT rat (rn7) fork — original 대비 수정 내역

upstream **copykat 1.1.0**, commit `b795ff793522499f814f6ae282aad1aab790902f`
(navinlabcode/copykat) 기준. 원본 함수 소스는 기존 Rocky9 CNV SIF에 설치된
패키지에서 `srcref`로 byte-faithful 복원했다(`docs/original_copykat_b795ff7.R`).
전체 통합 패치는 [`copykat_rn7.patch`](copykat_rn7.patch).

## 요약: 코드 3곳 + 데이터 1개

| # | 위치 | 변경 | 목적 |
| --- | --- | --- | --- |
| 1 | `copykat()` step 2 annotation 분기 | `else if(genome=="rn7")` 추가 | rat 값이 분기에 걸리게 (원래 미지원 값은 `else`가 없어 즉시 죽음) |
| 2 | `copykat()` 후반 gene-space 블록 | `if(genome=="mm10")` → `if(genome=="mm10" || genome=="rn7")` | rat이 mouse gene-space 다운스트림(220kb bin 없음)을 그대로 상속 |
| 3 | 파일 끝 | `annotateGenes.rn7()` 신규 추가 | `annotateGenes.mm10()` 복사본, 기본값 `full.anno=full.anno.rn7` |
| 4 | `R/sysdata.rda` | `full.anno.rn7` 추가 | mRatBN7.2 rat annotation table baking (기존 `full.anno`/`full.anno.mm10`/`DNA.hg20`/`cyclegenes`는 그대로 유지) |

메타데이터: `DESCRIPTION`의 Version `1.1.0` → `1.1.0.9001`, `NAMESPACE`에 `export(annotateGenes.rn7)` 1줄 추가.

## 수정하지 **않은** 것 (의도적)

- **hg20 전용 단계는 손대지 않음.** cell-cycle/HLA 유전자 제거(`if(genome=="hg20")`)와
  220kb genomic bin 변환(`convert.all.bins.hg20`)은 이미 `genome=="hg20"`로만 게이트되어
  있어, rn7은 자동으로 이 단계들을 건너뛴다.
- **`mgi_symbol` 컬럼명 유지.** `copykat()` 후반부가 `rownames(mat.adj) <- anno.mat2$mgi_symbol`로
  이 이름을 직접 참조하므로, rat symbol을 담되 컬럼명을 바꾸면 코드가 깨진다. 그래서
  `full.anno.rn7`도 6번째 컬럼명을 문자 그대로 `mgi_symbol`로 두고 rat gene symbol을 넣었다.
- **알고리즘 로직(필터·smoothing·baseline·MCMC segmentation)은 원본 그대로.**

## 상세 diff

### 편집 1 — annotation 분기 (원본 line 40 근처)
```diff
   } else if(genome=="mm10"){
   anno.mat <- annotateGenes.mm10(mat = rawmat, ID.type = id.type) #SYMBOL or ENSEMBLE
   dim(rawmat)
+  } else if(genome=="rn7"){
+  anno.mat <- annotateGenes.rn7(mat = rawmat, ID.type = id.type) #SYMBOL or ENSEMBLE (rat mRatBN7.2)
+  dim(rawmat)
   }
```
원본은 `hg20`/`mm10` 두 분기뿐이고 `else`가 없어, `genome="rn7"`을 그냥 넣으면
`anno.mat`이 생성되지 않아 다음 줄(`anno.mat[order(...)]`)에서 "object not found"로 죽는다.

### 편집 2 — gene-space 출력 블록 (원본 line 558 근처)
```diff
-  if(genome=="mm10") {
+  if(genome=="mm10" || genome=="rn7") {
     uber.mat.adj <- data.matrix(results.com)
```
이 블록이 mouse 방식(gene-space) 최종 예측·`CNA_results.txt`·heatmap을 담당한다.
hg20 전용 220kb bin 경로(`if(genome=="hg20") convert.all.bins.hg20(...)`)는 위에서 이미
분기되므로, rn7은 이 블록만 타면 gene-space로 끝난다.

### 편집 3 — `annotateGenes.rn7()` 신규 (파일 끝에 추가)
`annotateGenes.mm10()`과 본문이 동일하고 기본 인자만 다르다:
```diff
-annotateGenes.mm10 <-
-function(mat, ID.type="S", full.anno=full.anno.mm10){
+annotateGenes.rn7 <-
+function(mat, ID.type="S", full.anno=full.anno.rn7){
```
`id.type="S"`는 `full.anno$mgi_symbol`로, `"E"`는 `full.anno$ensembl_gene_id`로
count matrix와 교집합을 잡고, 중복 `mgi_symbol`을 제거한 뒤 annotation+expression을 합친다.
`mgi_symbol` 컬럼명을 유지했기 때문에 함수 본문은 한 글자도 바꿀 필요가 없다.

### 데이터 4 — `full.anno.rn7` (R/sysdata.rda)
`build_full_anno_rn7.R`가 `genes.gtf`(mRatBN7.2, GCA_015227675.2)의 `feature=="gene"`에서
생성. `full.anno.mm10`과 동일한 7컬럼 스키마:

| 컬럼 | 내용 |
| --- | --- |
| `abspos` | 실제 누적 좌표(염색체별 max end offset + start). 내장 `full.anno.mm10`의 상수-abspos 결함(unique 20/137,030)을 재현하지 않음 — rn7은 25,268 unique / 25,302 |
| `chromosome_name` | 정수 인코딩 (상염색체 1–20, X=21; Y/MT/scaffold 제외) |
| `start_position`, `end_position` | 유전자 좌표 |
| `ensembl_gene_id` | Ensembl gene ID (`id.type="E"` 매칭) |
| `mgi_symbol` | **rat** gene symbol (컬럼명은 유지) |
| `band` | `NA` (시각화용, 계산 미사용) |

## 원본과 diff 재생성
```bash
# 기존 SIF에서 pristine 원본을 복원해 두었음:
#   docs/original_copykat_b795ff7.R
# 통합 패치:
#   docs/copykat_rn7.patch   (= diff -u original fork)
diff -u docs/original_copykat_b795ff7.R copykat-rn7/R/copykat.R
```

## 검증 상태
- T1: `annotateGenes.rn7` 정의, `formals(copykat)$genome`, sysdata baking 확인.
- T2: 합성 rat 행렬로 `copykat(genome="rn7")` 완주(gene-space 출력, 리터럴 클래스).
- T3: 실데이터 310셀(HN00283643, IMPHH78w-20D)로 완주 — 아래 참조.
  결과는 `review_required` (DNA 근거 교차검증 전).
