#!/usr/bin/env python3
"""Example message designed with Slack's block kit builder:
{
	"blocks": [
		{
			"type": "section",
			"text": {
				"text": "*<https://google.com|LOOK Project> daily summary*",
				"type": "mrkdwn"
			}
		},
		{
			"type": "section",
			"text": {
				"text": "{n_images} image{plural_images} and {n_targets} target{plural_targets} processed in the past 24 hours.",
				"type": "mrkdwn"
			}
		},
		{
			"type": "section",
			"text": {
				"type": "mrkdwn",
				"text": "C/2022 E3   r=13.1 ± 0.03\n22P   g 11.14 ± 0.12"
			}
		}
	]
}
"""
import os
import json
import argparse
import requests
from urllib.parse import quote
import astropy.units as u
from astropy.io import ascii
from astropy.time import Time
from sbpy.data import natural_sort_key
from quick_look import setup_logger

parser = argparse.ArgumentParser(epilog='Requires environment variables '
                                 'LOOK_SLACK_USER_POST_URL, and LOOK_SLACK_NOTIFICATIONS_POST_URL')
# parser.add_argument('--date', help='generate report for this date')
parser.add_argument('-n', action='store_true', dest='debug',
                    help='no-op mode for testing')
parser.add_argument('--user', action='store_true',
                    help='post message to user for testing')
parser.add_argument('--interval', type=int, default=24,
                    help='time interval (hours) over which to check for observations')
parser.add_argument('--verbose', '-v', action='store_true',
                    help='increase command-line verbosity')
args = parser.parse_args()

logger = setup_logger('slack-notification',
                      'INFO' if args.verbose else 'WARNING')

base_url = 'https://quick-look.mkelley.dev'
user_post_url = os.getenv('LOOK_SLACK_USER_POST_URL')
notifications_post_url = os.getenv('LOOK_SLACK_NOTIFICATIONS_POST_URL')


def target_link(target):
    return f"<{base_url}/target/{quote(target.replace('/', '_'))}|{target}>"


def format_message(interval):
    t_start = (Time.now() - interval * u.hr)
    data = ascii.read('file-summary.txt')
    data = data[data['obs date'] > t_start.isot]
    data['basename'] = [os.path.basename(f[:-8]) for f in data['filename']]
    n = len(data)
    logger.info('%d files', n)

    if n == 0:
        return {
            "blocks": [
                {
                    "type": "section",
                    "text": {
                        "text": f"*<{base_url}/|Quick LOOK> daily summary*",
                        "type": "mrkdwn"
                    }
                },
                {
                    "type": "section",
                    "text": {
                        "text": f"No observations processed in the past {interval} hours.",
                        "type": "mrkdwn"
                    }
                },
            ]
        }

    skip = {row['file']: row['reason']
            for row in ascii.read('phot-skip.list')}
    data['skip'] = [skip.get(f, '') for f in data['basename']]
    logger.info('%d files skipped', sum(data['skip'] != ''))

    with open('stack-clusters.json', 'r') as inf:
        clusters = json.load(inf)

    data['stack prefix'] = [clusters.get(row['filename'], [''])[0][:-9]
                            for row in data]
    logger.info('%d files stacked', sum(data['stack prefix'] != ''))

    phot = ascii.read('phot.txt')
    phot = phot[phot['date'] > t_start.iso]
    targets_processed = set(phot['target'])

    binned = ascii.read('phot-binned.txt')
    binned = binned[binned['date'] > t_start.iso]
    targets_binned = set(binned['target'])
    binned = {row['stack prefix']: row for row in binned}

    n_images = len(data)
    plural_images = '' if n_images == 1 else 's'

    n_processed = len(targets_processed)
    plural_processed = '' if n_processed == 1 else 's'

    n_binned = len(targets_binned)
    # plural_binned = '' if n_binned == 1 else 's'

    message = {
        "blocks": [
            {
                "type": "section",
                "text": {
                    "text": f"*<{base_url}/|Quick LOOK> daily summary*",
                    "type": "mrkdwn"
                }
            },
            {
                "type": "section",
                "text": {
                    "text": f"{n_images} image{plural_images} and {n_processed} target{plural_processed} processed in the past {interval} hours.",
                    "type": "mrkdwn"
                }
            },
        ]
    }

    # summarize photometry
    rows = []
    for row in binned.values():
        rows.append(
            f"{target_link(row['target'])} — {row['catalog filter']} — {round(row['m5'], 2)} ± {round(row['merr5'], 2)}")

    if len(rows) > 0:
        text = "Photometry\n\n" + "\n".join(rows)
        if len(text) > 3000:
            text = text[:2997] + '...'
        message['blocks'].append({
            "type": "section",
            "text": {
                "text": text,
                "type": "mrkdwn"
            }
        })
    else:
        message['blocks'].append({
            "type": "section",
            "text": {
                "text": "No photometry",
                "type": "mrkdwn"
            }
        })

    # summarized processed but not binned
    rows = []
    for target in sorted(targets_processed - targets_binned, key=natural_sort_key):
        rows.append(target_link(target))
    if len(rows) > 0:
        text = "Processed, but failed photometry\n\n" + "\n".join(rows)
        if len(text) > 3000:
            text = text[:2997] + '...'
        message['blocks'].append({
            "type": "section",
            "text": {
                "text": text,
                "type": "mrkdwn"
            }
        })

    # summarized skipped data
    not_reduced = data[data['stack prefix'] == '']
    skipped = not_reduced[not_reduced['skip'] != '']
    if len(skipped) > 0:
        rows = []
        for target in sorted(set(skipped['target']), key=natural_sort_key):
            s = skipped[skipped['target'] == target]
            n = len(s)
            reasons = ", ".join(set(s['skip']))
            rows.append(f"{target_link(target)} ({n} files): {reasons}")
        message['blocks'].append({
            "type": "section",
            "text": {
                "text": "\n".join(rows),
                "type": "mrkdwn"
            }
        })

    # summarize other
    other = not_reduced[not_reduced['skip'] == '']
    if len(other) > 0:
        rows = []
        for target in sorted(set(other['target']), key=natural_sort_key):
            s = other[other['target'] == target]
            n = len(s)
            rows.append(
                f"{target_link(target)} ({n} files): status unknown, potentially mid-processing")
        message['blocks'].append({
            "type": "section",
            "text": {
                "text": "\n".join(rows),
                "type": "mrkdwn"
            }
        })

    return message


message = format_message(args.interval)
logger.info(message)
if args.user:
    r = requests.post(user_post_url, data=json.dumps(message))
else:
    r = requests.post(notifications_post_url, data=json.dumps(message))
