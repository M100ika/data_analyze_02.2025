#microbiome_analysis.py
"""
Microbiome Analysis Library
===========================

Данная библиотека содержит функции для анализа данных микробиома: расчёт альфа- и бета-разнообразия,
визуализация относительной абундантности таксонов (stacked plots, box-/violinplots), построение диаграмм Венна
и визуализацию соотношения Firmicutes/Bacteroidetes.

Перед использованием рекомендуется настроить глобальные переменные (например, пути для сохранения файлов,
флаги SAFE_DATA и SHOW, а также PALETTE) в соответствии с вашими нуждами.
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import kruskal
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova, anosim
from matplotlib_venn import venn3, venn3_circles
import plotly.express as px
import plotly.graph_objects as go
from matplotlib.patches import Ellipse

# ============================================================================
# Глобальные настройки (настройте под себя)
# ============================================================================

SAFE_DATA = False  # Если True, то результаты сохраняются в файлы
SHOW = False       # Если True, то графики отображаются

PATH_TO_ALPHA_OUTPUT_DATA = './alpha_output_data'
PATH_TO_ALPHA_OUTPUT_FIGURE = './alpha_output_figures'
PATH_TO_RELATIVE_ABUNDANCE_OUTPUT_DATA = './relative_abundance_output_data'
PATH_TO_RELATIVE_ABUNDANCE_OUTPUT_FIGURE = './relative_abundance_output_figures'
PATH_TO_VENN_OUTPUT_FIGURE = './venn_output_figures'
PATH_TO_BETA_OUTPUT_DATA = './beta_output_data'
PATH_TO_BETA_OUTPUT_FIGURE = './beta_output_figures'
PATH_TO_FB_OUTPUT_FIGURE = './fb_output_figures'
PATH_TO_FB_OUTPUT_DATA = './fb_output_data'
NAME_DF = 'default'
GROUP_LIST = 'default'
# PALETTE – словарь (например, {'Group1': 'red', 'Group2': 'blue', ...}). Задайте его самостоятельно.
PALETTE = {}

# ============================================================================
# Вспомогательные функции для обработки данных
# ============================================================================

def melt_to_pivot(df):
    """
    Преобразует DataFrame из длинного формата в сводную таблицу,
    где индексом являются ['#SampleID', 'GROUP'].
    """
    pivot_df = df.pivot(index=['#SampleID', 'GROUP'], columns='Taxonomy', values='RelativeAbundance').reset_index()
    pivot_df.fillna(0, inplace=True)
    pivot_df.columns.name = None
    return pivot_df

def filter_level(df, level):
    """Фильтрация данных по указанному уровню таксономии."""
    return df[df['Taxonomy'].str.contains(level)]

def trim_taxonomy_name(df):
    """Оставляет в колонке 'Taxonomy' только последний уровень таксономии (после разделения по '|')."""
    return df['Taxonomy'].str.split('|').str[-1]

def trim_taxonomy_name_2(df, col_idx=0):
    """Оставляет в указанной колонке только последний уровень таксономии."""
    return df.iloc[:, col_idx].str.split('|').str[-1]

def trim_phylum(df):
    """Удаляет префикс 'p__' и оставляет только название филума."""
    df['Taxonomy'] = df['Taxonomy'].str.replace('p__', '').str.split('|').str[0]
    return df

def save_results_to_excel(results, file_path):
    """
    Сохраняет результаты ANOVA и теста Тьюки в Excel-файл с разными листами.
    """
    with pd.ExcelWriter(file_path, engine='openpyxl') as writer:
        for key, data in results.items():
            if 'ANOVA' in data and data['ANOVA'] is not None:
                data['ANOVA'].to_excel(writer, sheet_name=f'Group_{key}_ANOVA')
            
            if 'Tukey' in data and data['Tukey'] is not None:
                data['Tukey'].to_excel(writer, sheet_name=f'Group_{key}_Tukey')

# ============================================================================
# Функции для расчёта альфа-разнообразия
# ============================================================================

def alpha_calculation(pivot_df, normalize=True):
    """
    Вычисляет индексы альфа-разнообразия (shannon, simpson, observed, pielou_e)
    на основе сводной таблицы pivot_df, 
    где строки = сэмплы, столбцы = таксоны + два служебных столбца: '#SampleID', 'GROUP'.

    Параметры
    ---------
    pivot_df : pd.DataFrame
        Сводная таблица, полученная из melt_to_pivot. 
        Индексы/строки = сэмплы, столбцы = таксоны (и '#SampleID','GROUP').
    normalize : bool
        Если True, тогда каждую строку нормируем, чтобы сумма значений по таксонам = 1.
        Это нужно, если `pivot_df` содержит сырые (ненормированные) счёты.

    Возвращает
    ---------
    alpha_diversity_df : pd.DataFrame
        Таблица со столбцами: ['SampleID', 'GROUP', 'shannon', 'simpson', 'observed', 'pielou_e'].
    """

    meta_cols = ['#SampleID', 'GROUP']
    # Вырезаем таксоны, оставляем только счётчики (или уже относительные).
    pivot_data = pivot_df.drop(columns=meta_cols, errors='ignore').fillna(0)

    if normalize:
        # Нормируем каждую строчку, чтобы сумма была 1 (если она > 0).
        row_sums = pivot_data.sum(axis=1)
        pivot_data = pivot_data.div(row_sums, axis=0).fillna(0)

    data = pivot_data.values

    # Собираем идентификаторы
    ids = pivot_df['#SampleID'].values
    groups = pivot_df['GROUP'].values

    # Считаем альфа-метрики с помощью scikit-bio
    shannon_diversity = alpha_diversity('shannon', data, ids)
    simpson_diversity = alpha_diversity('simpson', data, ids)
    observed_otus = alpha_diversity('observed_otus', data, ids)
    pielou_e = alpha_diversity('pielou_e', data, ids)

    alpha_diversity_df = pd.DataFrame({
        'SampleID': ids,
        'GROUP': groups,
        'shannon': shannon_diversity,
        'simpson': simpson_diversity,
        'observed': observed_otus,
        'pielou_e': pielou_e
    })
    alpha_diversity_df = alpha_diversity_df.reset_index(drop=True).fillna(0)

    # Сохраняем, если нужно
    if SAFE_DATA:
        alpha_diversity_df.to_excel(f'{PATH_TO_ALPHA_OUTPUT_DATA}/{NAME_DF}_alpha_general.xlsx', index=False)
    return alpha_diversity_df

def compare_alpha_diversity_groups(alpha_diversity_df, pairs_dict):
    """
    Сравнивает значения индекса shannon между группами (по словарю пар).
    
    Параметры:
      - alpha_diversity_df: DataFrame с рассчитанными индексами альфа-разнообразия.
      - pairs_dict: словарь, где ключ – идентификатор, а значение – список групп для сравнения.
    
    Возвращает словарь с результатами ANOVA и, при наличии значимой разницы, результатами теста Тьюки.
    """
    results = {}
    
    for key, groups in pairs_dict.items():
        group_data = alpha_diversity_df[alpha_diversity_df['GROUP'].isin(groups)]
        
        model = ols('shannon ~ C(GROUP)', data=group_data).fit()
        anova_result = sm.stats.anova_lm(model, typ=2)
        
        if anova_result['PR(>F)'][0] < 0.05:
            tukey_result = pairwise_tukeyhsd(endog=group_data['shannon'], groups=group_data['GROUP'], alpha=0.05)
            tukey_summary = pd.DataFrame(data=tukey_result._results_table.data[1:], 
                                         columns=tukey_result._results_table.data[0])
            results[key] = {'ANOVA': anova_result, 'Tukey': tukey_summary}
        else:
            results[key] = {'ANOVA': anova_result, 'Tukey': None}

    return results


def alpha_stats_tukey(alpha_diversity_df, index, subgroup_name):
    """
    Выполняет анализ дисперсии (ANOVA) и тест Тьюки для указанного индекса альфа-разнообразия.

    Возвращает: (anova_table, tukey_results_df).
    """

    # Модель ANOVA
    model = ols(f'{index} ~ C(GROUP)', data=alpha_diversity_df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)

    # Пост-хок Тьюки (учитывает множественные сравнения между группами)
    tukey = pairwise_tukeyhsd(endog=alpha_diversity_df[index],
                              groups=alpha_diversity_df['GROUP'], alpha=0.05)
    tukey_results_df = pd.DataFrame(data=tukey._results_table.data[1:], 
                                    columns=tukey._results_table.data[0])
    # При желании можно сохранять
    if SAFE_DATA:
        tukey_results_df.to_excel(f'{PATH_TO_ALPHA_OUTPUT_DATA}/{subgroup_name}_{index}_alpha_stats_tukey.xlsx', index=False)

    return anova_table, tukey_results_df


def p_value_alpha(alpha_diversity_df, index):
    """
    Проводит Тьюки (post-hoc) для пар групп по указанному индексу (index).
    Возвращает DataFrame (матрица p-value) размером G x G, где G — количество групп.

    Внимание:
    pairwise_tukeyhsd сам по себе уже учитывает множественные сравнения между группами.
    """
    tukey_result = pairwise_tukeyhsd(
        endog=alpha_diversity_df[index],
        groups=alpha_diversity_df['GROUP'],
        alpha=0.05
    )

    # Переводим выход Тьюки в DataFrame
    tukey_summary_df = pd.DataFrame(data=tukey_result._results_table.data[1:], 
                                    columns=tukey_result._results_table.data[0])

    # Создаём квадратную матрицу p-value
    groups = alpha_diversity_df['GROUP'].unique()
    p_value_matrix = pd.DataFrame(index=groups, columns=groups, data=1.0)

    for _, row in tukey_summary_df.iterrows():
        g1 = row['group1']
        g2 = row['group2']
        p_val = row['p-adj']
        p_value_matrix.loc[g1, g2] = p_val
        p_value_matrix.loc[g2, g1] = p_val

    if SAFE_DATA:
        p_value_matrix.to_excel(f'{PATH_TO_ALPHA_OUTPUT_DATA}/{NAME_DF}_{index}_pvals_matrix.xlsx', index=True)
    return p_value_matrix


def alpha_boxplot(alpha_diversity_df, p_value_matrix, alpha_index, subgroup_name, group_order, figsize=(4, 5)):
    """
    Строит boxplot для индекса альфа-разнообразия alpha_index с аннотацией значимости.
    """
    # Проверка на пустые данные
    if alpha_diversity_df.empty:
        print(f"Warning: No data for {alpha_index} in the provided DataFrame.")
        return

    plt.figure(figsize=figsize)

    # Проверка, что все группы присутствуют в данных
    existing_groups = alpha_diversity_df['GROUP'].unique()
    group_order = [group for group in group_order if group in existing_groups]

    if not group_order:
        print(f"Warning: No valid groups found for {alpha_index}.")
        return

    # 1) Рисуем boxplot
    ax = sns.boxplot(
        x='GROUP',
        y=alpha_index,
        data=alpha_diversity_df,
        order=group_order,
        boxprops=dict(facecolor='white', edgecolor='black', linewidth=2),
        whiskerprops=dict(color='black', linewidth=2),
        capprops=dict(color='black', linewidth=2),
        medianprops=dict(color='black', linewidth=2),
        flierprops=dict(marker=''),
        width=0.6
    )

    # 2) Накладываем stripplot/swarmplot вручную
    for group in group_order:
        sub_df = alpha_diversity_df[alpha_diversity_df['GROUP'] == group]
        if not sub_df.empty:
            sns.stripplot(
                x='GROUP', y=alpha_index,
                data=sub_df,
                order=[group],
                color=PALETTE.get(group, 'gray'),
                marker='o', size=5,
                jitter=True,
                ax=ax
            )

    plt.xticks(rotation=90)

    plt.title(f'{alpha_index.capitalize()} Diversity Index by Group')
    plt.xlabel('Group')
    plt.ylabel(alpha_index.capitalize())

    # 3) Аннотация p-values
    group_combinations = [
        (group_order[i], group_order[j])
        for i in range(len(group_order))
        for j in range(i+1, len(group_order))
    ]
    y_max = alpha_diversity_df[alpha_index].max()
    y_min = alpha_diversity_df[alpha_index].min()
    y_range = y_max - y_min
    y_offset = y_max + 0.05 * y_range
    y_step = 0.05 * y_range

    annotation_count = 0
    for g1, g2 in group_combinations:
        p_val = p_value_matrix.loc[g1, g2]
        if p_val < 0.05:
            x1 = group_order.index(g1)
            x2 = group_order.index(g2)
            y = y_offset + annotation_count * y_step
            plt.plot([x1, x1, x2, x2],
                     [y, y + y_step/2, y + y_step/2, y],
                     lw=1.5, c='black')
            p_val_text = f"p={p_val:.3e}" if p_val >= 1e-3 else "p<0.001"
            plt.text((x1 + x2) / 2, y + y_step/2, p_val_text,
                     ha='center', va='bottom', color='black')
            annotation_count += 1

    plt.tight_layout()
    if SAFE_DATA:
        plt.savefig(f'{PATH_TO_ALPHA_OUTPUT_FIGURE}/{subgroup_name}_{alpha_index}_boxplot.pdf', bbox_inches='tight')
    if SHOW:
        plt.show()
    plt.close()


def alpha_all(alpha_df, subgroup_name, figsize=[(4,5),(4,5),(4,5),(4,5)]):
    """
    Запуск всех этапов для 4 метрик: 
    - вычисление матрицы p-value (Тьюки) 
    - ANOVA + Тьюки (в логи)
    - построение boxplot'ов
    """
    for metric, size in zip(['shannon','simpson','observed','pielou_e'], figsize):
        # 1) Получаем матрицу p-value
        pval_mat = p_value_alpha(alpha_df, metric)
        
        # 2) Выполняем ANOVA + Tukey (просто чтобы сохранить результат)
        anova_tab, tukey_df = alpha_stats_tukey(alpha_df, metric, subgroup_name)
        print(f"\n=== {metric.upper()} ANOVA ===\n", anova_tab)
        print(f"=== {metric.upper()} Tukey ===\n", tukey_df.head())

        # 3) Строим boxplot
        alpha_boxplot(alpha_df, pval_mat, metric, subgroup_name, GROUP_LIST, figsize=size)
        

def alpha_violinplot(alpha_diversity_df, p_value, alpha_index, group_order):
    """
    Строит violin plot для выбранного индекса альфа-разнообразия с аннотацией значимости.
    
    Параметры аналогичны функции alpha_boxplot.
    """
    plt.figure(figsize=(8, 6))
    markers = ['o', 's', '^', 'D'] 

    sns.violinplot(x='GROUP', y=alpha_index, data=alpha_diversity_df,
                   order=group_order,
                   palette=PALETTE,
                   cut=0,
                   scale='width',
                   inner=None)

    sns.stripplot(x='GROUP', y=alpha_index, data=alpha_diversity_df,
                  order=group_order,
                  hue='GROUP',
                  palette=PALETTE,
                  dodge=False,
                  marker='o',
                  size=6,
                  jitter=True)

    plt.title(f'{alpha_index.capitalize()} Diversity Index by Group')
    plt.xlabel('Group')
    plt.ylabel(alpha_index.capitalize())

    groups = group_order
    group_combinations = [(groups[i], groups[j]) for i in range(len(groups))
                          for j in range(i+1, len(groups))]

    y_max = alpha_diversity_df[alpha_index].max()
    y_range = y_max - alpha_diversity_df[alpha_index].min()
    y_offset = y_max + y_range * 0.05
    y_step = y_range * 0.05

    for idx, (group1, group2) in enumerate(group_combinations):
        p_val = p_value.loc[group1, group2]
        x1 = group_order.index(group1)
        x2 = group_order.index(group2)
        y = y_offset + idx * y_step
        plt.plot([x1, x1, x2, x2], [y, y + y_step/2, y + y_step/2, y],
                 lw=1.5, c='black')
        p_val_text = f"p = {p_val:.4f}" if p_val >= 0.001 else "p < 0.001"
        plt.text((x1 + x2) / 2, y + y_step/2 + y_range*0.01,
                 p_val_text, ha='center', va='bottom', color='black')

    plt.legend([], [], frameon=False)
    plt.tight_layout()
    if SAFE_DATA:
        plt.savefig(f'{PATH_TO_ALPHA_OUTPUT_FIGURE}/{NAME_DF}_{alpha_index}_violinplot.pdf',
                    bbox_inches='tight')
    if SHOW:
        plt.show()
    plt.close()

# ============================================================================
# Функции для построения stacked plots и других графиков таксономической абундантности
# ============================================================================

def stack_plot_test(taxon_abundance, pairs, name_suffix='', orientation='vertical'):
    """Строит stacked bar plot для топ-10 таксонов с использованием matplotlib.
    Параметры:
    - taxon_abundance: DataFrame с колонками 'GROUP', 'Taxonomy' и 'Value'.
    - pairs: список групп для включения.
    - name_suffix: суффикс для имен файлов.
    - orientation: 'vertical' или 'horizontal'.
    Возвращает словарь с p-value для каждого таксона."""
    
    # Словарь для замены имен групп
    group_rename_dict = {
        'CHJ-BPD/DS+': 'VFGM-BPD/DS',
        'CHJ-BPD/DS-': 'OB/SD'
    }
    
    # Отладка: проверяем уникальные группы
    print("Исходные уникальные группы:", taxon_abundance['GROUP'].unique())
    print("Группы в pairs:", pairs)
    
    # Создаем копию данных и очищаем от лишних пробелов
    taxon_abundance_copy = taxon_abundance.copy()
    taxon_abundance_copy['GROUP'] = taxon_abundance_copy['GROUP'].str.strip()
    
    # Заменяем имена групп
    taxon_abundance_copy['GROUP'] = taxon_abundance_copy['GROUP'].replace(group_rename_dict)
    
    # Отладка: проверяем группы после замены
    print("Группы после замены:", taxon_abundance_copy['GROUP'].unique())
    
    # Также обновляем список pairs с новыми именами (и очищаем от пробелов)
    cleaned_pairs = [str(group).strip() for group in pairs]
    updated_pairs = [group_rename_dict.get(group, group) for group in cleaned_pairs]
    
    print("Обновленные pairs:", updated_pairs)
    
    filtered_abundance = taxon_abundance_copy[taxon_abundance_copy['GROUP'].isin(updated_pairs)]
    
    # Остальной код остается без изменений...
    mean_abundance = filtered_abundance.groupby(['GROUP', 'Taxonomy'])['Value'].mean().reset_index()
    mean_abundance['Value'] = pd.to_numeric(mean_abundance['Value'], errors='coerce')
    top_taxa = mean_abundance.groupby('Taxonomy')['Value'].sum().nlargest(10).index
    top_taxon_abundance = mean_abundance[mean_abundance['Taxonomy'].isin(top_taxa)]
    
    p_values = {}
    for taxon in top_taxa:
        groups = filtered_abundance[filtered_abundance['Taxonomy'] == taxon].groupby('GROUP')['Value'].apply(list)
        if len(groups) > 1:
            _, p_value = kruskal(*groups)
            p_values[taxon] = p_value
        else:
            p_values[taxon] = None
    
    legend_labels = [f"{taxon} {'*' * (3 if p_values[taxon] is not None and p_values[taxon] < 0.001 else 2 if p_values[taxon] is not None and p_values[taxon] < 0.01 else 1 if p_values[taxon] is not None and p_values[taxon] < 0.05 else 0)}"
                     for taxon in top_taxa]
    
    pivot_table = top_taxon_abundance.pivot_table(index='GROUP', columns='Taxonomy',
                                                  values='Value', aggfunc='mean')
    pivot_table = pivot_table.reindex(updated_pairs)
    sorted_taxa = list(top_taxa)
    pivot_table = pivot_table[sorted_taxa]
    
    colors = sns.color_palette('deep', n_colors=len(top_taxa))
    taxa_colors = dict(zip(top_taxa, colors))
    
    fig, ax = plt.subplots(figsize=(2, 5))
    plot_kind = 'barh' if orientation == 'horizontal' else 'bar'
    pivot_table.loc[:, top_taxa].plot(kind=plot_kind, stacked=True,
                                      color=[taxa_colors[taxon] for taxon in top_taxa],
                                      width=1, edgecolor='white', ax=ax)
    
    ax.set_title(f'Top 10 Taxa {name_suffix}')
    if orientation == 'vertical':
        ax.set_ylabel('Percentage of Relative Abundance')
        ax.set_xlabel('Groups')
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    else:
        ax.set_xlabel('Percentage of Relative Abundance')
        ax.set_ylabel('Groups')
        ax.set_yticklabels(ax.get_yticklabels(), rotation=90)
    
    ax.set_facecolor('none')
    fig.set_facecolor('none')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    patches = [plt.Rectangle((0,0),1,1, color=taxa_colors[taxon]) for taxon in top_taxa]
    plt.legend(patches, legend_labels, title="Taxa Categories",
               bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)
    
    plt.tight_layout()
    
    if SAFE_DATA:
        pivot_table.to_excel(f'{PATH_TO_RELATIVE_ABUNDANCE_OUTPUT_DATA}/{NAME_DF}_{name_suffix}_relative_abundance.xlsx', index=True)
        plt.savefig(f'{PATH_TO_RELATIVE_ABUNDANCE_OUTPUT_FIGURE}/{NAME_DF}_stack_{name_suffix}.pdf', bbox_inches='tight')
        with open(f'{PATH_TO_RELATIVE_ABUNDANCE_OUTPUT_DATA}/p_values_{name_suffix}.txt', 'w') as file:
            file.write(f'p_values_{name_suffix}: \n{p_values}\n\n')
    
    if SHOW:
        plt.show()
    plt.close()
    
    return p_values


def stack_plot_test_plotly(taxon_abundance, pairs, name_suffix='', orientation='vertical'):
    """
    Строит stacked bar plot для топ-10 таксонов с использованием Plotly.
    
    Параметры аналогичны функции stack_plot_test.
    Возвращает словарь с p-value для каждого таксона.
    """
    filtered_abundance = taxon_abundance[taxon_abundance['GROUP'].isin(pairs)]
    mean_abundance = filtered_abundance.groupby(['GROUP', 'Taxonomy'])['Value'].mean().reset_index()
    mean_abundance['Value'] = pd.to_numeric(mean_abundance['Value'], errors='coerce')

    top_taxa = mean_abundance.groupby('Taxonomy')['Value'].sum().nlargest(50).index
    top_taxon_abundance = mean_abundance[mean_abundance['Taxonomy'].isin(top_taxa)]

    p_values = {}
    for taxon in top_taxa:
        groups = filtered_abundance[filtered_abundance['Taxonomy'] == taxon].groupby('GROUP')['Value'].apply(list)
        if len(groups) > 1:
            _, p_value = kruskal(*groups)
            p_values[taxon] = p_value
        else:
            p_values[taxon] = None

    legend_labels = [f"{taxon} {'*' * (3 if p_values[taxon] is not None and p_values[taxon] < 0.001 else 2 if p_values[taxon] is not None and p_values[taxon] < 0.01 else 1 if p_values[taxon] is not None and p_values[taxon] < 0.05 else 0)}"
                     for taxon in top_taxa]

    pivot_table = top_taxon_abundance.pivot_table(index='GROUP', columns='Taxonomy',
                                                  values='Value', aggfunc='mean')
    pivot_table = pivot_table[top_taxa]
    pivot_table = pivot_table.reindex(pairs)

    if orientation == 'vertical':
        fig = go.Figure()
        for taxon in top_taxa:
            fig.add_trace(go.Bar(
                x=pivot_table.index,
                y=pivot_table[taxon],
                name=f"{taxon} {'*' * (3 if p_values[taxon] is not None and p_values[taxon] < 0.001 else 2 if p_values[taxon] is not None and p_values[taxon] < 0.01 else 1 if p_values[taxon] is not None and p_values[taxon] < 0.05 else 0)}",
                marker=dict(color=px.colors.qualitative.Plotly[top_taxa.tolist().index(taxon) % len(px.colors.qualitative.Plotly)])
            ))
        fig.update_layout(barmode='stack', title=f'Top 10 Taxa {name_suffix}',
                          yaxis_title='Percentage of Relative Abundance', xaxis_title='Groups')
    else:
        fig = go.Figure()
        for taxon in top_taxa:
            fig.add_trace(go.Bar(
                y=pivot_table.index,
                x=pivot_table[taxon],
                name=f"{taxon} {'*' * (3 if p_values[taxon] is not None and p_values[taxon] < 0.001 else 2 if p_values[taxon] is not None and p_values[taxon] < 0.01 else 1 if p_values[taxon] is not None and p_values[taxon] < 0.05 else 0)}",
                orientation='h',
                marker=dict(color=px.colors.qualitative.Plotly[top_taxa.tolist().index(taxon) % len(px.colors.qualitative.Plotly)])
            ))
        fig.update_layout(barmode='stack', title=f'Top 10 Taxa {name_suffix}',
                          xaxis_title='Percentage of Relative Abundance', yaxis_title='Groups')

    fig.update_layout(
        barmode='stack',
        title=f'Top 10 Taxa {name_suffix}',
        xaxis_title='Percentage of Relative Abundance' if orientation == 'horizontal' else 'Groups',
        yaxis_title='Groups' if orientation == 'horizontal' else 'Percentage of Relative Abundance',
        xaxis=dict(showgrid=False, gridcolor='LightGray', gridwidth=0.5),
        yaxis=dict(showgrid=False, gridcolor='LightGray', gridwidth=0.5),
        plot_bgcolor='white',
        legend_title="Taxa Categories",
        width=600,
        height=800
    )

    if SAFE_DATA:
        pivot_table.to_excel(f'{PATH_TO_RELATIVE_ABUNDANCE_OUTPUT_DATA}/{NAME_DF}_{name_suffix}_plotly.xlsx', index=True)
        fig.write_image(f'{PATH_TO_RELATIVE_ABUNDANCE_OUTPUT_FIGURE}/{NAME_DF}_plotly_{name_suffix}.pdf')
        fig.write_html(f'{PATH_TO_RELATIVE_ABUNDANCE_OUTPUT_FIGURE}/{NAME_DF}_plotly_{name_suffix}.html')
        with open(f'{PATH_TO_RELATIVE_ABUNDANCE_OUTPUT_DATA}/p_values_plotly_{name_suffix}.txt', 'w') as file:
            file.write(f'p_values_{name_suffix}: \n{p_values}\n\n')

    if SHOW:
        fig.show()

    return p_values

def mean_plotbars(df, name_suffix, group_orders, relative='RelativeAbundance', height=0.1):
    """
    Строит столбчатую диаграмму (bar plot) с error bars и аннотациями значимости для топ-10 таксонов.
    
    Параметры:
      - df: DataFrame с колонками 'Taxonomy', 'GROUP' и relative (например, 'RelativeAbundance').
      - name_suffix: суффикс для сохранения графика.
      - group_orders: порядок групп для отображения.
      - relative: имя колонки с относительной абундантностью.
      - height: смещение для аннотаций.
    
    Возвращает DataFrame с объединёнными стандартными отклонениями.
    """
    std_devs_per_taxonomy = df.groupby('Taxonomy')[relative].std().fillna(0)
    df_with_std_dev = df.merge(std_devs_per_taxonomy.rename('StdDev'), on='Taxonomy', how='left')
    top_taxa = df.groupby('Taxonomy')[relative].mean().nlargest(10).index
    df_with_std_dev = df_with_std_dev[df_with_std_dev['Taxonomy'].isin(top_taxa)]

    plt.figure(figsize=(10, 9))
    sns.barplot(
        x='Taxonomy',
        y=relative,
        hue='GROUP',
        data=df_with_std_dev,
        palette=PALETTE,
        alpha=0.5,
        capsize=0.1,
        errwidth=1,
        width=0.8,
        hue_order=group_orders
    )

    for i, taxonomy in enumerate(df_with_std_dev['Taxonomy'].unique()):
        bar_pos = np.array(range(len(df_with_std_dev['Taxonomy'].unique())))
        plt.errorbar(
            x=bar_pos[i],
            y=df_with_std_dev[df_with_std_dev['Taxonomy'] == taxonomy][relative].mean(),
            yerr=std_devs_per_taxonomy[taxonomy],
            fmt='none',
            c='k',
            capsize=3,
        )
        background_color = 'white' if i % 2 else 'gray'
        plt.axvspan(bar_pos[i] - 0.5, bar_pos[i] + 0.5, color=background_color, alpha=0.1)

        data = [df_with_std_dev.loc[(df_with_std_dev['Taxonomy'] == taxonomy) &
                                    (df_with_std_dev['GROUP'] == group), relative].values 
                for group in df_with_std_dev['GROUP'].unique()]
        _, p_value = kruskal(*data)
        
        if relative == 'RelativeAbundance':
            plt.text(bar_pos[i], df_with_std_dev[relative].max() - height,
                     f"{'*' * (3 if p_value < 0.001 else 2 if p_value < 0.01 else 1 if p_value < 0.05 else 0)}",
                     ha='center', va='bottom', color='black', fontsize=10)
        else:
            plt.text(bar_pos[i], df_with_std_dev[relative].max() - 1.08,
                     f"{'*' * (3 if p_value < 0.001 else 2 if p_value < 0.01 else 1 if p_value < 0.05 else 0)}",
                     ha='center', va='bottom', color='black', fontsize=10)

    custom_dash = [20, 5, 20, 5]
    plt.axhline(0, color='black', linewidth=1, dashes=custom_dash, alpha=0.5)

    plt.xticks(rotation=60, ha='right')
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.xlabel('Taxonomy')
    plt.ylabel('Relative Abundance')
    plt.title('Relative Abundance of Taxonomies in Different Groups')
    plt.legend(title='GROUP', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(color='gray', linestyle='--', linewidth=0.1)
    plt.tight_layout()
    if SAFE_DATA:
        plt.savefig(f'{PATH_TO_RELATIVE_ABUNDANCE_OUTPUT_FIGURE}/mean_plot_{name_suffix}_{relative}.pdf', bbox_inches='tight')
    if SHOW:
        plt.show()
    plt.close()
    return df_with_std_dev

def venn3_plot(taxa_1, taxa_2, taxa_3, group_tuple):
    """
    Строит диаграмму Венна для трёх наборов таксонов.
    
    Параметры:
      - taxa_1, taxa_2, taxa_3: наборы (или списки) таксонов.
      - group_tuple: кортеж с названиями групп.
    """
    plt.figure(figsize=(8, 8))
    v = venn3([set(taxa_1), set(taxa_2), set(taxa_3)], group_tuple)
    plt.title("Venn Diagram of Taxonomies in Different Groups")

    segments = ['100', '010', '001', '110', '101', '011', '111']
    colors = ['skyblue', 'lightgreen', 'salmon', 'blue', 'green', 'red', 'purple']
    for segment, color in zip(segments, colors):
        patch = v.get_patch_by_id(segment)
        if patch is not None:
            patch.set_color(color)
            patch.set_alpha(0.4)

    venn3_circles(subsets=(1, 1, 1, 1, 1, 1, 1), linestyle='dashed')
    for label in v.subset_labels:
        if label:
            label.set_fontsize(14)
            label.set_fontweight('bold')
    for label, text in zip(v.set_labels, group_tuple):
        label.set_fontsize(16)
        label.set_fontweight('bold')
        label.set_color('red')

    if SAFE_DATA:
        plt.savefig(f'{PATH_TO_VENN_OUTPUT_FIGURE}/venn_{GROUP_LIST}.pdf', bbox_inches='tight')
    if SHOW:
        plt.show()
    plt.close()

# ============================================================================
# Функция для визуализации бета-разнообразия
# ============================================================================

def beta_diversity_all(long_format_df,  subgroup_name = ''):
    beta_diversity_visualization(long_format_df, beta_diversity_method="braycurtis", subgroup_name = subgroup_name)
    beta_diversity_visualization(long_format_df, beta_diversity_method="jaccard",  subgroup_name = subgroup_name)
    beta_diversity_visualization(long_format_df, beta_diversity_method='euclidean',  subgroup_name = subgroup_name)
    beta_diversity_visualization(long_format_df, beta_diversity_method='canberra',  subgroup_name = subgroup_name)
    beta_diversity_visualization(long_format_df, beta_diversity_method='cosine',  subgroup_name = subgroup_name)



def beta_diversity_visualization(long_format_df, beta_diversity_method="braycurtis", subgroup_name = ''):
    """
    Визуализирует бета-разнообразие с помощью PCoA, отображая группы с эллипсами.
    
    Параметры:
      - long_format_df: DataFrame в длинном формате с колонками '#SampleID', 'Taxonomy',
        'RelativeAbundance' и 'GROUP'.
      - beta_diversity_method: метод расчёта расстояний (по умолчанию "braycurtis").
    """
    data_pivot = long_format_df.pivot_table(index='#SampleID', columns='Taxonomy',
                                            values='RelativeAbundance', fill_value=0)
    samples = data_pivot.index.tolist()
    group_labels = long_format_df.drop_duplicates(subset=['#SampleID']).set_index('#SampleID')['GROUP']

    dm = beta_diversity(beta_diversity_method, data_pivot.values, ids=samples)
    dm_df = pd.DataFrame(dm.data, index=samples, columns=samples)
    
    np.nan_to_num(dm.data, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    pcoa_results = pcoa(dm)
    pcoa_df = pcoa_results.samples
    pcoa_df['GROUP'] = group_labels.values

    def draw_ellipse(position, covariance, ax=None, **kwargs):
        ax = ax or plt.gca()
        if covariance.shape == (2, 2):
            U, s, Vt = np.linalg.svd(covariance)
            angle = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
            width, height = 2 * np.sqrt(s)
        else:
            angle = 0
            width = height = 2 * np.sqrt(covariance)
        for nsig in range(1, 2):
            ell = Ellipse(xy=position, width=nsig * width, height=nsig * height,
                          angle=angle, **kwargs, linewidth=3)
            ax.add_patch(ell)

    permanova_results = permanova(dm, group_labels)
    anosim_results = anosim(dm, group_labels)
        
    plt.figure(figsize=(9, 6))
    sns.scatterplot(x='PC1', y='PC2', hue='GROUP', data=pcoa_df, s=50, alpha=0.75,
                    palette=PALETTE, edgecolor="black", linewidth=0.5)
    for label, group_df in pcoa_df.groupby("GROUP"):
        color = PALETTE.get(label, 'blue')
        draw_ellipse(group_df[["PC1", "PC2"]].mean(), group_df[["PC1", "PC2"]].cov(),
                     alpha=0.25, facecolor='None', edgecolor=color)

    anosim_text = f"ANOSIM R = {anosim_results['test statistic']:.3f}, p-value = {anosim_results['p-value']:.3f}"
    permanova_text = f"PERMANOVA R = {permanova_results['test statistic']:.3f}, p-value = {permanova_results['p-value']:.3f}"
    plt.text(.94, 0.032, anosim_text, ha="right", va="bottom",
             transform=plt.gca().transAxes, fontsize=12)
    plt.text(.94, 0.002, permanova_text, ha="right", va="bottom",
             transform=plt.gca().transAxes, fontsize=12)
    plt.title(f"{subgroup_name} - PCoA plot with {beta_diversity_method} dissimilarity", fontsize=18)
    plt.xlabel(f"PC1 ({pcoa_results.proportion_explained.iloc[0]:.2%} explained)", fontsize=14)
    plt.ylabel(f"PC2 ({pcoa_results.proportion_explained.iloc[1]:.2%} explained)", fontsize=14)
    plt.axhline(0, color='grey', linestyle='--')
    plt.axvline(0, color='grey', linestyle='--')
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.legend(title='GROUP', fontsize=12, title_fontsize=12)
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    if SAFE_DATA:
        with open(f'{PATH_TO_BETA_OUTPUT_DATA}/{NAME_DF}_{beta_diversity_method}_stats.txt', 'w') as file:
            file.write(f'{beta_diversity_method} PERMANOVA results: \n{permanova_results}\n\n'
                       f'{beta_diversity_method} ANOSIM results: \n{anosim_results}')
        dm_df.to_excel(f'{PATH_TO_BETA_OUTPUT_DATA}/{subgroup_name}_{beta_diversity_method}.xlsx', index=True)
        plt.savefig(f'{PATH_TO_BETA_OUTPUT_FIGURE}/{subgroup_name}_{beta_diversity_method}.pdf', bbox_inches='tight')
    if SHOW:
        plt.show()
    plt.close()

# ============================================================================
# Функция для визуализации соотношения Firmicutes/Bacteroidetes
# ============================================================================

def fb_ratio_viz(df, group_order, firmicutes = 'p__Firmicutes', bacteroidetes = 'p__Bacteroidota'):
    """
    Визуализирует соотношение F/B (Firmicutes/Bacteroidetes) по группам.
    
    Параметры:
      - df: DataFrame с колонками '#SampleID', 'GROUP', 'Taxonomy' и 'RelativeAbundance'.
      - group_order: список групп в нужном порядке.
    """
    from scipy.stats import kruskal  # убедимся, что функция kruskal импортирована

    def calculate_fb_ratio(df):
        df_firmicutes = df[df['Taxonomy'].str.contains(firmicutes)] 
        df_bacteroidetes = df[df['Taxonomy'].str.contains(bacteroidetes)]
        if len(df_firmicutes) == 0 or len(df_bacteroidetes) == 0:
            print('Impossible to create figure')
            print('len df_firmicutes:', len(df_firmicutes))
            print('len df_bacteroidetes:', len(df_bacteroidetes))
        firmicutes_sum = df_firmicutes.groupby(['#SampleID', 'GROUP'])['RelativeAbundance'].sum().reset_index()
        bacteroidetes_sum = df_bacteroidetes.groupby(['#SampleID', 'GROUP'])['RelativeAbundance'].sum().reset_index()
        fb_combined = pd.merge(firmicutes_sum, bacteroidetes_sum, on=['#SampleID', 'GROUP'],
                               suffixes=('_firmicutes', '_bacteroidetes'))
        fb_combined['F/B Ratio'] = fb_combined['RelativeAbundance_firmicutes'] / fb_combined['RelativeAbundance_bacteroidetes']
        return fb_combined

    def calculate_p_values(fb_data):
        groups = fb_data['GROUP'].unique()
        p_values = {}
        for i in range(len(groups)):
            for j in range(i+1, len(groups)):
                group1 = groups[i]
                group2 = groups[j]
                _, p_val = kruskal(
                    fb_data[fb_data['GROUP'] == group1]['F/B Ratio'],
                    fb_data[fb_data['GROUP'] == group2]['F/B Ratio']
                )
                p_values[(group1, group2)] = p_val
        return p_values

    fb_combined = calculate_fb_ratio(df)
    
    # Фильтрация выбросов (пример для определённых групп)
    fb_combined = fb_combined[~((fb_combined['GROUP'] == 'CDAI — II') & (fb_combined['F/B Ratio'] > 6000))]
    fb_combined = fb_combined[~((fb_combined['GROUP'] == 'CDAI — I') & (fb_combined['F/B Ratio'] > 180))]
    
    p_values = calculate_p_values(fb_combined)

    plt.figure(figsize=(6, 6))
    ax = sns.boxplot(data=fb_combined, x='GROUP', y='F/B Ratio',
                     palette=['white']*len(group_order), order=group_order,
                     saturation=1, width=0.4)
    markers = ['o', 's', '^'] 
    for i in range(3):
        sns.swarmplot(x='GROUP', y='F/B Ratio', data=fb_combined,
                      hue='GROUP', palette=PALETTE, marker=markers[i], size=6)
    
    plt.title('F/B Ratio Across Groups')
    y_max = ax.get_ylim()[1] * 1.05

    for (group1, group2), p_val in p_values.items():
        x1 = group_order.index(group1)
        x2 = group_order.index(group2)
        ax.plot([x1, x2], [y_max, y_max], color='black')
        ax.text((x1 + x2) / 2, y_max, f'{p_val:.3f}', ha='center', va='bottom')
        y_max *= 1.05  

    plt.ylabel('F/B Ratio')
    
    if SAFE_DATA:
        plt.savefig(f'{PATH_TO_FB_OUTPUT_FIGURE}/fb_ratio_{NAME_DF}.pdf', bbox_inches='tight')
    if SHOW:
        plt.show()
    plt.close()

    fb_combined['log(F/B Ratio)'] = np.log(fb_combined['F/B Ratio'] + 1e-10)
    fb_combined = fb_combined[np.isfinite(fb_combined['log(F/B Ratio)'])]

    plt.figure(figsize=(10, 6))
    for group in group_order:
        sns.kdeplot(fb_combined[fb_combined['GROUP'] == group]['log(F/B Ratio)'],
                    fill=True, label=group, color=PALETTE.get(group, None))
    plt.title('Density of log(Firmicutes/Bacteroidetes Ratio) Across Groups')
    plt.xlabel('log(Firmicutes/Bacteroidetes Ratio)')
    plt.ylabel('Density')
    plt.legend(title='Group')
    if SAFE_DATA:
        fb_combined.to_excel(f'{PATH_TO_FB_OUTPUT_DATA}/fb_combined_{NAME_DF}.xlsx', index=True)
        plt.savefig(f'{PATH_TO_FB_OUTPUT_FIGURE}/Density_{NAME_DF}.pdf', bbox_inches='tight')
    if SHOW:
        plt.show()
    plt.close()

# ============================================================================
# Конец библиотеки
# ============================================================================

if __name__ == '__main__':
    # Пример использования:
    # Здесь можно разместить тестовые вызовы функций при запуске этого файла напрямую.
    print("Модуль microbiome_analysis импортирован. Импортируйте функции для дальнейшего анализа.")
