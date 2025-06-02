#lib.py
import numpy as np
import pandas as pd

GROUP = None 


def update_asv_columns(merged_otu_df, biom_data):
    """
    Функция заменяет названия столбцов ASV в merged_otu_df на таксономию из biom_data.
    
    Аргументы:
      merged_otu_df: DataFrame, в котором первые столбцы (например, "#SampleID", "GROUP")
                     остаются без изменений, а остальные (ASV1, ASV2, …) будут заменены.
      biom_data: Словарь, загруженный из OTU.biom (формат JSON), где под ключом 'rows'
                 находится список записей с ключами 'id' и 'metadata'.
                 
    Возвращает:
      DataFrame с обновлёнными названиями столбцов.
    """

    asv_to_tax = {}
    for entry in biom_data.get('rows', []):
        asv_id = entry.get('id')
        taxonomy_levels = entry.get('metadata', {}).get('taxonomy', [])
        taxonomy_str = '|'.join(taxonomy_levels)
        asv_to_tax[asv_id] = taxonomy_str


    new_columns = []
    for col in merged_otu_df.columns:

        if col in asv_to_tax:
            new_columns.append(asv_to_tax[col])
        else:
            new_columns.append(col)
            
    merged_otu_df.columns = new_columns
    return merged_otu_df


def read_data(file_path, sep='\t'):
    """Функция для чтения данных из файла."""
    return pd.read_csv(file_path, sep=sep)

def preprocess_meta_map(meta_map):
    """Преобразование данных meta_map."""
    meta_map = meta_map.drop(columns=['SequencingRun'])
    for column in meta_map.select_dtypes(include=['object']):  
        meta_map[column] = meta_map[column].str.replace('.extendedFrags.fastq.gz', '', regex=False)
        #meta_map[column] = meta_map[column].str.replace('^Oral', '', regex=True)
    #meta_map['fastqFile'] = meta_map['fastqFile'].apply(lambda x: x.split('_')[0])
    meta_map = meta_map.rename(columns={'fastqFile':'ID'})
    #meta_map['ID'] = pd.to_numeric(meta_map['ID'], errors='coerce')
    return meta_map

def merge_and_prepare_datasets(meta_map, meta, df, replace_dict = None):
    """Слияние и подготовка наборов данных."""
    global GROUP
    merged_df = pd.merge(meta_map, meta, on='ID')
    final_meta_df = merged_df[['#SampleID', 'ID', GROUP]]
    final_meta_df = final_meta_df.rename(columns={GROUP:'GROUP'})
    final_meta_df = final_meta_df[final_meta_df['GROUP'].notna() & (final_meta_df['GROUP'] != '')]

    df_transposed = df.T.reset_index()
    df_transposed.columns = ['#SampleID'] + list(df.iloc[:,0])
    df_transposed = df_transposed.drop(index=0)
    
    merged_df = pd.merge(final_meta_df, df_transposed, on='#SampleID', how='left').drop(columns=['ID'])
    
    if replace_dict != None:
        merged_df['GROUP'] = merged_df['GROUP'].replace(replace_dict)
    
    return final_meta_df, merged_df

def melt_and_prepare_taxonomy(merged_df):
    """Преобразование формата данных и подготовка таксономии."""
    # Переводим таблицу в длинный формат
    long_format_df = pd.melt(merged_df, id_vars=['#SampleID', 'GROUP'], 
                               var_name='Taxonomy', value_name='Value')
    
    # Добавляем столбец OTU с формированием идентификатора
    long_format_df['OTU'] = 'ASV' + (long_format_df.index + 1).astype(str).str.zfill(4)
    long_format_df = long_format_df[['OTU', 'Taxonomy', '#SampleID', 'GROUP', 'Value']]
    
    # Применяем функцию transform_taxonomy, затем удаляем лишний префикс "k_"
    long_format_df['Taxonomy'] = (long_format_df['Taxonomy']
                                  .apply(transform_taxonomy)
                                  .str.replace(r'^k_', '', regex=True)  # удаляем ведущий k_
                                  .replace(r'k_noHit\|p_', 'noHit', regex=True))
    
    # Применяем дополнительное преобразование таксономии, если необходимо
    long_format_df['Taxonomy'] = long_format_df['Taxonomy'].apply(trim_taxonomy)
    
    # Фильтруем записи, где таксономия не равна "noHit"
    return long_format_df[long_format_df['Taxonomy'] != 'noHit']


def transform_taxonomy(row):
    """Трансформация таксономии."""
    tax_levels = row.split(';')
    prefixes = ['k_', 'p_', 'c_', 'o_', 'f_', 'g_', 's_']
    transformed = '|'.join([f"{prefix}{level}" for prefix, level in zip(prefixes, tax_levels)])
    return transformed

def trim_taxonomy(taxonomy):
    """Обрезка таксономии."""
    pos = taxonomy.find('_?')
    return taxonomy[:taxonomy.rfind('|', 0, pos)] if pos != -1 else taxonomy

def calculate_relative_abundance(df):
    """Расчёт относительного изобилия таксонов."""
    return rel_abundance(df)

def filter_level(df, level):
    """Фильтрация данных до уровня рода."""
    return df[df['Taxonomy'].str.contains(level)]

def rel_abundance(df):
    """Расчет относительного изобилия."""
    taxon_abundance = df.groupby(['#SampleID', 'Taxonomy', 'GROUP'])['Value'].sum().reset_index()
    total_abundance = taxon_abundance.groupby(['#SampleID', 'GROUP'])['Value'].sum().reset_index().rename(columns={'Value': 'Total'})
    taxon_abundance = pd.merge(taxon_abundance, total_abundance, on=['#SampleID', 'GROUP'])
    taxon_abundance['RelativeAbundance'] = taxon_abundance.apply(lambda row: row['Value'] / row['Total'] if row['Total'] > 0 else 0, axis=1)
    return taxon_abundance

def remove_low_variance_taxa(df, threshold=0.01):
    """
    Фильтрация таксонов на основе пороговой вариативности относительного изобилия.
    
    Параметры:
    df (pd.DataFrame): Датафрейм с относительным изобилием таксонов.
    threshold (float): Пороговая вариативность для фильтрации таксонов.
    
    Возвращает:
    pd.DataFrame: Отфильтрованный датафрейм.
    """

    variance = df.groupby('Taxonomy')['RelativeAbundance'].var()

    taxa_above_threshold = variance[variance > threshold].index

    filtered_df = df[df['Taxonomy'].isin(taxa_above_threshold)]
    
    return filtered_df


def filter_threshold(df, threshold=0.01):
    """
    Фильтрация таксонов на основе порогового значения относительного изобилия.
    
    Параметры:
    df (pd.DataFrame): Датафрейм с относительным изобилием таксонов.
    threshold (float): Пороговое значение для фильтрации таксонов.
    
    Возвращает:
    pd.DataFrame: Отфильтрованный датафрейм.
    """

    sum_abundance = df.groupby('Taxonomy')['RelativeAbundance'].sum()
    
    taxa_above_threshold = sum_abundance[sum_abundance > threshold].index
    
    filtered_df = df[df['Taxonomy'].isin(taxa_above_threshold)]
    
    return filtered_df

def safe_data(df, path_out, gr_name, safe_file = True):
    pivot_table_lefse = df.pivot_table(
        index='Taxonomy', columns=['GROUP', '#SampleID'],
        values='RelativeAbundance', aggfunc='sum'
    ).fillna(0)
    pivot_table_lefse = pivot_table_lefse.reset_index().T.reset_index().T

    if safe_file:
        pivot_table_lefse.rename(columns={'Taxonomy':'GROUP'}).to_csv(
            f'{path_out}/lefse/data/{gr_name}.tsv',
            sep='\t', index=False, header=False
        )
    return pivot_table_lefse


def final_lefse_df(df, path_out, gr_name, safe_file = False):
    pivot_table_lefse = df.pivot_table(index='Taxonomy', columns=['GROUP', '#SampleID'], values='RelativeAbundance', aggfunc='sum').fillna(0)
    pivot_table_lefse = pivot_table_lefse.reset_index().T.reset_index().T
    
    pivot_table_lefse = safe_data(df, path_out, gr_name, safe_file = safe_file)
   
    df_genus = filter_level(df, 'g_')
    safe_data(df_genus, path_out, gr_name+'_genus', safe_file = safe_file)
    df_species = filter_level(df, 's_')
    safe_data(df_species, path_out, gr_name+'_species', safe_file = safe_file)
    return pivot_table_lefse

def clean_df(df):
        # Remove columns where names end with "?"
        df = df.loc[:, ~df.columns.str.endswith('?')]
        # Shorten column names by keeping only the last part after the last ";"
        df.columns = [name.split(';')[-1] if ';' in name else name for name in df.columns]
        df = df.drop(columns = ['GROUP',''])
        
        return df


def to_maaslin_data(df, meta_df):
    
    df_transposed = df.T.reset_index()
    df_transposed.columns = ['#SampleID'] + list(df['Genus'])
    df_transposed = df_transposed.drop(df_transposed.index[0])
    merged_df = pd.merge(meta_df, df_transposed, on='#SampleID', how='left')
    merged_df = merged_df.drop(columns=['ID'])
    to_maaslin = clean_df(merged_df)
    return to_maaslin


def attach_subgroup_column(safe_df, meta_df, meta_id='ID', df_id='#SampleID', subgroup_col='Subgroup'):
    """
    Добавляет столбец Subgroup (или другой, если задать subgroup_col) в датафрейм safe_df,
    используя данные из meta_df.
    
    Параметры:
    ----------
    safe_df : pd.DataFrame
        Таблица в "длинном" формате с колонкой df_id (обычно '#SampleID').
    meta_df : pd.DataFrame
        Таблица с колонками [meta_id, subgroup_col], где meta_id (обычно 'ID') – 
        идентификатор, соответствующий #SampleID, а subgroup_col – нужная колонка.
    meta_id : str
        Название столбца в meta_df, соответствующее #SampleID в safe_df.
    df_id : str
        Название столбца в safe_df, по которому будем объединять (обычно '#SampleID').
    subgroup_col : str
        Название столбца в meta_df, который хотим добавить (по умолчанию 'Subgroup').
    
    Возвращает:
    ----------
    pd.DataFrame
        Копию safe_df, дополненную колонкой subgroup_col.
    """

    # Создадим копию meta_df и переименуем 'ID' -> '#SampleID', чтобы они совпали
    meta_renamed = meta_df[[meta_id, subgroup_col]].rename(columns={meta_id: df_id})

    # Мержим safe_df с meta_renamed
    merged_df = pd.merge(safe_df, meta_renamed, on=df_id, how='left')
    
    # Перемещаем колонку 'Subgroup' сразу после 'GROUP'
    cols = merged_df.columns.tolist()
    group_idx = cols.index('GROUP')  # Находим индекс столбца 'GROUP'
    cols.insert(group_idx + 1, cols.pop(cols.index(subgroup_col)))  # Перемещаем 'Subgroup'
    merged_df = merged_df[cols]

    return merged_df

