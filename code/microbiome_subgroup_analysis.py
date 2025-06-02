"""
microbiome_subgroup_analysis.py

Библиотека для анализа микробиома с учётом 'Subgroup'.
-----------------------------------------------------------------
Она не меняет функций из microbiome_analysis.py,
а лишь добавляет обёртки (wrappers), которые:

1) Фильтруют DataFrame по нужному Subgroup.
2) Вызывают функции из microbiome_analysis.py
   (например, alpha_all, stack_plot_test, mean_plotbars).
3) Повторяют процедуру для нужного списка Subgroup.

Пример использования:

import microbiome_subgroup_analysis as msa

# Предположим, у вас есть DF с колонками:
# ['#SampleID','Taxonomy','GROUP','Subgroup','RelativeAbundance', ...]

# Хотим проанализировать только 'MAM-J':
alpha_df = msa.analyze_subgroup_alpha(DF, 'MAM-J', group_list=['OB-BPD/DS ','VFGM-BPD/DS ','VFGM','CN/SD '])

# Если нужно запустить сразу весь пакет анализов (alpha, stacked, mean bars):
alpha_df, pvals_stack, df_meanbars = msa.analyze_subgroup_all(DF, 'MAM-J', group_list=GROUP_LIST)

-----------------------------------------------------------------
"""

import pandas as pd
import numpy as np
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

import microbiome_analysis as ma

def analyze_subgroup_alpha(df, subgroup_name, group_list=None):
    """
    Проводит альфа-диверсификационный анализ (alpha_all) только по подмножеству, 
    где df['Subgroup'] == subgroup_name.

    Параметры
    ---------
    df : pd.DataFrame
        Должен содержать колонки:
        - '#SampleID' (строка)
        - 'GROUP' (строка)
        - 'Subgroup' (строка)
        - 'Taxonomy', 'RelativeAbundance' (или другое имя, но по умолчанию 'RelativeAbundance')
    subgroup_name : str
        Название подгруппы, например: 'MAM-J', 'MAM-I' и т.д.
    group_list : list or None
        Список групп, которые нужно оставить. Если None, берутся все группы
        (какие есть в подмножестве).

    Возвращает
    ---------
    alpha_df : pd.DataFrame
        Результаты alpha_calculation (Shannon, Simpson, Observed, Pielou_e) 
        по подмножеству с нужным Subgroup.
    """
    # 1. Фильтруем по Subgroup
    df_sub = df[df['Subgroup'] == subgroup_name].copy()
    if group_list is not None:
        df_sub = df_sub[df_sub['GROUP'].isin(group_list)]

    # 2. Преобразуем в pivot (по #SampleID, GROUP)
    pivot_df = ma.melt_to_pivot(df_sub)
    # 3. Считаем альфа-диверсити
    alpha_df = ma.alpha_calculation(pivot_df, normalize=True)
    # 4. Вызываем ma.alpha_all, чтобы автоматически построить boxplot и т.д.
    #    (если SHOW=True, то отобразит, если SAFE_DATA=True — сохранит)
    ma.alpha_all(alpha_df, subgroup_name)

    return alpha_df


def analyze_subgroup_stackplot(df, subgroup_name, group_list=None, name_suffix='stack_subgroup'):
    """
    Строит stacked bar plot (stack_plot_test) для выбранного Subgroup.
    Возвращает p-values (Kruskal), которые считает stack_plot_test.

    Параметры
    ---------
    df : pd.DataFrame
        Должен содержать 'GROUP','Subgroup','Taxonomy','RelativeAbundance' как минимум.
    subgroup_name : str
        Какой Subgroup брать (например, 'MAM-J').
    group_list : list or None
        Какие GROUP использовать. Если None, берём все, какие в подмножестве.
    name_suffix : str
        Суффикс для имени файла, если SAFE_DATA=True (stack_subgroup по умолчанию).

    Возвращает
    ---------
    p_values : dict
        Словарь {taxon: p_value (Kruskal)}
    """
    df_sub = df[df['Subgroup'] == subgroup_name].copy()
    if group_list is not None:
        df_sub = df_sub[df_sub['GROUP'].isin(group_list)]

    # stack_plot_test ожидает DataFrame c колонками ['GROUP','Taxonomy','Value']
    df_stack = df_sub[['GROUP','Taxonomy','RelativeAbundance']].copy()
    df_stack.rename(columns={'RelativeAbundance': 'Value'}, inplace=True)

    if group_list is None:
        # Берём уникальные группы, которые остались
        pairs = sorted(df_stack['GROUP'].unique().tolist())
    else:
        pairs = group_list

    # Вызываем функцию из microbiome_analysis
    p_values = ma.stack_plot_test(
        taxon_abundance=df_stack,
        pairs=pairs,
        name_suffix=f"{name_suffix}_{subgroup_name}",
        orientation='vertical'
    )
    return p_values


def analyze_subgroup_mean_bars(df, subgroup_name, group_list=None, name_suffix='mean_subgroup'):
    """
    Строит столбчатый график (mean_plotbars) для top-10 таксонов, только по строкам с заданным Subgroup.

    Параметры
    ---------
    df : pd.DataFrame
        Должен содержать 'GROUP','Subgroup','Taxonomy','RelativeAbundance' 
    subgroup_name : str
        Какой Subgroup берем
    group_list : list or None
        Какие GROUP оставляем
    name_suffix : str
        Суффикс для имени файла

    Возвращает
    ---------
    df_meanbars : pd.DataFrame
        То, что возвращает mean_plotbars (таблица с добавленным столбцом 'StdDev'),
        но только по выбранному Subgroup.
    """
    df_sub = df[df['Subgroup'] == subgroup_name].copy()
    if group_list is not None:
        df_sub = df_sub[df_sub['GROUP'].isin(group_list)]

    if group_list is None:
        group_orders = sorted(df_sub['GROUP'].unique().tolist())
    else:
        group_orders = group_list

    df_meanbars = ma.mean_plotbars(
        df=df_sub,
        name_suffix=f"{name_suffix}_{subgroup_name}",
        group_orders=group_orders,
        relative='RelativeAbundance',
        height=0.1
    )
    return df_meanbars


def analyze_subgroup_all(df, subgroup_name, group_list=None):
    """
    Запускает основные этапы сразу:
    1) alpha-анализ (alpha_all)
    2) stack plot
    3) mean bars
    4) beta
    5) Kruskal–Wallis по каждому таксону

    и возвращает результаты.
    """
    # 1) alpha
    alpha_df = analyze_subgroup_alpha(df, subgroup_name, group_list=group_list)
    
    # 2) stack
    pvals_stack = analyze_subgroup_stackplot(df, subgroup_name, group_list=group_list)
    
    # 3) mean bars
    df_mean = analyze_subgroup_mean_bars(df, subgroup_name, group_list=group_list)

    # 4) beta
    ma.beta_diversity_all(df[df['Subgroup'] == subgroup_name], subgroup_name)

    # 5) Kruskal–Wallis по каждому таксону
    df_kruskal = analyze_subgroup_kruskal_significant_taxa(df,
                                                           subgroup_name=subgroup_name,
                                                           group_list=group_list)
    
    return alpha_df, pvals_stack, df_mean, df_kruskal


def analyze_all_subgroups(df, subgroup_list, group_list=None):
    """
    Запускает analyze_subgroup_all(...) по очереди для каждого Subgroup из списка.
    
    Возвращает:
      results : dict
         ключ = название Subgroup,
         значение = tuple (alpha_df, pvals_stack, df_mean, df_kruskal)
    """
    results = {}
    for sg in subgroup_list:
        print(f"\n=== Analyzing Subgroup: {sg} ===")
        alpha_df, pvals_stack, df_mean, df_kruskal = analyze_subgroup_all(
            df, sg, group_list=group_list
        )
        results[sg] = (alpha_df, pvals_stack, df_mean, df_kruskal)
    return results


def test_taxa_kruskal_fdr(df, group_col='GROUP', value_col='RelativeAbundance', taxonomy_col='Taxonomy'):
    """
    Выполняет Kruskal–Wallis тест по каждому таксону (таксон = groupby по taxonomy_col),
    смотрит различия между группами (group_col). Корректирует p-value методом FDR.
    Возвращает DataFrame:
        [Taxonomy, p_value_raw, p_value_fdr],
    отсортированный по p_value_fdr.
    """
    results = []
    unique_taxa = df[taxonomy_col].unique()

    # Узнаём список групп
    groups_in_data = df[group_col].unique()

    for taxon in unique_taxa:
        df_taxon = df[df[taxonomy_col] == taxon]
        grouped_data = [
            df_taxon.loc[df_taxon[group_col] == g, value_col].dropna().values
            for g in groups_in_data
        ]
        # Проверим, что в каждой группе есть хоть какие-то значения
        if all(len(arr) > 0 for arr in grouped_data):
            stat, p_val = kruskal(*grouped_data)
            results.append((taxon, p_val))

    if len(results) > 0:
        taxa, pvals = zip(*results)
        pvals_corrected = multipletests(pvals, alpha=0.05, method='fdr_bh')[1]
        results_df = pd.DataFrame({
            taxonomy_col: taxa,
            'p_value_raw': pvals,
            'p_value_fdr': pvals_corrected
        })
        return results_df.sort_values('p_value_fdr')
    else:
        return pd.DataFrame(columns=[taxonomy_col, 'p_value_raw', 'p_value_fdr'])


def analyze_subgroup_kruskal_significant_taxa(df, subgroup_name, 
                                              group_list=None,
                                              abundance_col='RelativeAbundance',
                                              taxonomy_col='Taxonomy'):
    """
    Для указанного subgroup_name:
      1) Фильтрует df по Subgroup == subgroup_name.
      2) (Опционально) фильтрует по списку групп group_list.
      3) Запускает Kruskal–Wallis тест через test_taxa_kruskal_fdr, 
         чтобы узнать, какие таксоны (taxonomy_col) статистически отличаются между группами.
      4) Возвращает DataFrame с колонками [Taxonomy, p_value_raw, p_value_fdr],
         отсортированный по p_value_fdr.
    
    Если SAFE_DATA == True, то результат сохранится в Excel/CSV.
    
    Пример использования:
        df_kruskal = analyze_subgroup_kruskal_significant_taxa(DF, 'MAM-J')
        df_kruskal_signif = df_kruskal[df_kruskal['p_value_fdr'] < 0.05]
    """
    df_sub = df[df['Subgroup'] == subgroup_name].copy()
    if group_list is not None:
        df_sub = df_sub[df_sub['GROUP'].isin(group_list)]
    
    if df_sub.empty:
        print(f"[WARNING] Нет данных для Subgroup={subgroup_name}. Пропуск Kruskal-теста.")
        return pd.DataFrame(columns=[taxonomy_col, 'p_value_raw', 'p_value_fdr'])
    
    # Запуск теста
    df_kruskal = test_taxa_kruskal_fdr(
        df_sub,
        group_col='GROUP',
        value_col=abundance_col,
        taxonomy_col=taxonomy_col
    )
    
    # Сохранение, если нужно
    from microbiome_analysis import SAFE_DATA, PATH_TO_RELATIVE_ABUNDANCE_OUTPUT_DATA
    if SAFE_DATA and not df_kruskal.empty:
        out_path = f"{PATH_TO_RELATIVE_ABUNDANCE_OUTPUT_DATA}/{subgroup_name}_kruskal_significant_taxa.xlsx"
        df_kruskal.to_excel(out_path, index=False)
        print(f"[INFO] Kruskal results for '{subgroup_name}' saved to:", out_path)
    
    return df_kruskal
