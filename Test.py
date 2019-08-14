import scholarly
from typing import Dict, List
import os
import sys

import http.client
import urllib
from bs4 import BeautifulSoup
import re

import pybtex


class GoogleScholarSearch:
    """
    @brief This class searches Google Scholar (http://scholar.google.com)

    Search for articles and publications containing terms of interest.

    Usage example:\n
    <tt>
    > from google_search import *\n
    > searcher = GoogleScholarSearch()\n
    > searcher.search(['breast cancer', 'gene'])
    </tt>
    """

    def __init__(self):
        """
        @brief Empty constructor.
        """
        self.SEARCH_HOST = "scholar.google.com"
        self.SEARCH_BASE_URL = "/scholar"

    def search(self, terms, limit=10):
        """
        @brief This function searches Google Scholar using the specified terms.

        Returns a list of dictionarys. Each
        dictionary contains the information related to the article:
            "URL"		: link to the article/n
            "Title"		: title of the publication/n
            "Authors"	: authors (example: DF Easton, DT Bishop, D Ford)/n
            "JournalYear" 	: journal name & year (example: Nature, 2001)/n
            "JournalURL"	: link to the journal main website (example: www.nature.com)/n
            "Abstract"	: abstract of the publication/n
            "NumCited"	: number of times the publication is cited/n
            "Terms"		: list of search terms used in the query/n

        @param terms List of search terms
        @param limit Maximum number of results to be returned (default=10)
        @return List of results, this is the empty list if nothing is found
        """
        params = urllib.parse.urlencode({'q': "+".join(terms), 'num': limit})
        headers = {'User-Agent': 'Mozilla/4.0 (compatible; MSIE 5.5; Windows NT)'}

        url = self.SEARCH_BASE_URL + "?" + params
        conn = http.client.HTTPConnection(self.SEARCH_HOST)
        conn.request("GET", url, {}, headers)
        print(url)
        resp = conn.getresponse()

        if resp.status == 200:
            html = resp.read()
            results = []
            html = html.decode('ascii', 'ignore')


            # Screen-scrape the result to obtain the publication information
            soup = BeautifulSoup(html, 'html.parser')
            citations = 0
            for record in soup('p', {'class': 'g'}):
                # Includeds error checking
                topPart = record.first('span', {'class': 'w'})

                pubURL = topPart.a['href']
                # Clean up the URL, make sure it does not contain '\' but '/' instead
                pubURL = pubURL.replace('\\', '/')

                pubTitle = ""

                for part in topPart.a.contents:
                    pubTitle += str(part)

                if pubTitle == "":
                    match1 = re.findall('<b>\[CITATION\]</b></font>(.*)- <a', str(record))
                    match2 = re.split('- <a', match1[citations])
                    pubTitle = re.sub('</?(\S)+>', "", match2[0])
                    citations = citations + 1

                authorPart = record.first('font', {'color': 'green'}).string
                if str(authorPart) == 'Null':
                    authorPart = ''
                    # Sometimes even BeautifulSoup can fail, fall back to regex
                    m = re.findall('<font color="green">(.*)</font>', str(record))
                    if len(m) > 0:
                        authorPart = m[0]
                num = authorPart.count(" - ")
                # Assume that the fields are delimited by ' - ', the first entry will be the
                # list of authors, the last entry is the journal URL, anything in between
                # should be the journal year
                idx_start = authorPart.find(' - ')
                idx_end = authorPart.rfind(' - ')
                pubAuthors = authorPart[:idx_start]
                pubJournalYear = authorPart[idx_start + 3:idx_end]
                pubJournalURL = authorPart[idx_end + 3:]
                # If (only one ' - ' is found) and (the end bit contains '\d\d\d\d')
                # then the last bit is journal year instead of journal URL
                if pubJournalYear == '' and re.search('\d\d\d\d', pubJournalURL) is not None:
                    pubJournalYear = pubJournalURL
                    pubJournalURL = ''

                # This can potentially fail if all of the abstract can be contained in the space
                # provided such that no '...' is found
                delimiter = soup.firstText("...").parent
                pubAbstract = ""
                while str(delimiter) != 'Null' and (str(delimiter) != '<b>...</b>' or pubAbstract == ""):
                    pubAbstract += str(delimiter)
                    delimiter = delimiter.nextSibling
                pubAbstract += '<b>...</b>'

                match = re.search('Cited by ([^<]*)', str(record))
                pubCitation = ''
                if match != None:
                    pubCitation = match.group(1)
                results.append({
                    "URL": pubURL,
                    "Title": pubTitle,
                    "Authors": pubAuthors,
                    "JournalYear": pubJournalYear,
                    "JournalURL": pubJournalURL,
                    "Abstract": pubAbstract,
                    "NumCited": pubCitation,
                    "Terms": terms
                })
            return results
        else:
            print("ERROR: ")
            print(resp.status, resp.reason)
            return []


def read_listing(listing_path: str) -> Dict[int, str]:
    result = {}
    with open(listing_path, 'r', encoding="utf8") as listing_file:
        for line in listing_file:
            clean_line = line.strip()
            if clean_line != '':
                items = clean_line.split('\t', 1)
                result[int(items[0].strip('[').strip(']'))] = items[1]
    return result


def test_scholarly(item: str) -> List[Dict]:
    results = []
    for item in scholarly.search_pubs_query(item):
        results.append(item.bib)
    return results


def test_google_search(item: str) -> List[Dict]:
    results = []
    search_item = GoogleScholarSearch()
    for item in search_item.search(item):
        results.append(item)
    return results


def attempt_internet_pull(folder_path: str):
    listing_path = os.path.join(folder_path, "Bib_159.txt")
    out_listing = os.path.join(folder_path, "Bib_159.bib")
    bib_items = read_listing(listing_path)
    with open(out_listing, 'w', encoding="utf8") as write_file:
        for bib_number, bib_text in bib_items.items():
            scholar_bibs = test_scholarly(bib_text)
            #scholar_bibs = test_google_search([bib_text])
            print("{}) {} - results: {}".format(bib_number, len(scholar_bibs), bib_text))
            if len(scholar_bibs) != 1:
                print("item {} has {} scholar listings".format(bib_number, len(scholar_bibs)), file=sys.stderr)
            for res in scholar_bibs:
                print(res)
                bib_dict = res.copy()
                bib_dict['index'] = bib_number
                write_file.write('{}\n'.format(bib_dict))
            print("\n")


def print_bib_IEEE(folder_path: str):
    listing_path = os.path.join(folder_path, "Full_bibliography.bib")
    #pybtex.format_from_file(listing_path, .style.pybtexplain)


def set_couter_tags(file_path: str):
    IMPORTANT_RE = r"@(?P<type>.+){(?P<name>.+),"
    counter = 160
    with open(file_path,'r', encoding="utf8") as input_file, open("{}_counter".format(file_path),'w', encoding="utf8") as output_file:
        for line in input_file:
            matcher = re.search(IMPORTANT_RE, line)
            if matcher:
                type = matcher.group("type")
                output_file.write("@{}{}{:03d},\n".format(type, '{', counter))
                counter -= 1
            else:
                output_file.write(line)


if __name__ == "__main__":
    FOLDER = "D:\\Matan\\Dropbox\\PHd\\Final Thesis\\data\\"
    #attempt_internet_pull(FOLDER)
    #print_bib_IEEE(FOLDER)
    set_couter_tags(os.path.join(FOLDER, "Full_bibliography.bib"))
